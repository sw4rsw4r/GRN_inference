from networkx import *
import torch
import torch_geometric.transforms as T
from torch_geometric.nn.conv import GATv2Conv
from torch_geometric.data import Data
from torch_geometric.utils import negative_sampling
from torch_geometric.utils import to_networkx
from torch_geometric.datasets import Planetoid
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
import pickle


def load_embedding(path_files, TEST_ID, type):
    dt_embeddings = pd.read_csv(path_files+'/' + TEST_ID + '_embeddings_for_'+type+'.txt',
                                header=None, sep='\t')
    embeddings = np.asarray([x.astype(np.float32)
                            for x in dt_embeddings.values])
    return (embeddings)


def load_edge_index(path_files, TEST_ID, type):
    dt_edge_index = pd.read_csv(path_files+'/' + TEST_ID + '_edge_index_for_'+type+'.txt',
                                header=None, sep='\t')
    edge_index = torch.from_numpy(np.array(dt_edge_index).T)
    return (edge_index)


def build_graph(path_files, TEST_ID, type):
    edge_index = load_edge_index(path_files, TEST_ID, type)
    graph = Data(edge_index=edge_index.cuda())
    embeddings = load_embedding(path_files, TEST_ID, type)
    graph.num_nodes = len(embeddings)
    graph.x = torch.from_numpy(embeddings).type(torch.float32).cuda()
    return (graph)


class Net(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = GATv2Conv(in_channels, hidden_channels * 2, heads=2)
        self.Linear = torch.nn.Linear(hidden_channels * 4, hidden_channels * 4)
        self.conv2 = GATv2Conv(hidden_channels * 4, out_channels, heads=2)

    def encode(self, x, edge_index):
        x1 = self.conv1(x, edge_index).relu()
        x2 = self.Linear(x1)
        return self.conv2(x2, edge_index)

    def decode(self, z, edge_label_index):
        return (z[edge_label_index[0]] * z[edge_label_index[1]]).sum(
            dim=-1
        )  # product of a pair of nodes on each edge

    def decode_all(self, z):
        prob_adj = z @ z.t()
        return (prob_adj > 0).nonzero(as_tuple=False).t()


def train_link_predictor(
    model, train_data, val_data, test_data, optimizer, criterion, TEST_ID, rng_seed=1234, neg_ratio=1.0, file_suffix=''
):
    fname_output = 'GAT/data/'+TEST_ID+'_rng'+str(rng_seed)+'_neg'+str(neg_ratio) + file_suffix
    # Define the early stopping criteria
    best_val_auc = 0.
    validation_auc_list = []
    patience = 500.
    num_epoche_no_improve = 0
    with open( fname_output + '.log', 'w') as f:
        f.write('Start training..\n')
    n_epochs = 5000
    for epoch in range(1, n_epochs + 1):
        # set the model to training mode
        model.train()
        # reset the optimizer gradients
        optimizer.zero_grad()
        # encode the input data
        z = model.encode(train_data.x, train_data.edge_index)
        # sampling training negatives for every training epoch
        neg_edge_index = negative_sampling(
            edge_index=train_data.edge_index, num_nodes=train_data.num_nodes,
            num_neg_samples=train_data.edge_label_index.size(1), method='sparse')
        # concatenate the edge labels for positive and negative edges
        edge_label_index = torch.cat(
            [train_data.edge_label_index, neg_edge_index],
            dim=-1,
        )
        edge_label = torch.cat([
            train_data.edge_label,
            train_data.edge_label.new_zeros(neg_edge_index.size(1))
        ], dim=0)
        # decode the edge label indices
        out = model.decode(z, edge_label_index).view(-1)
        # compute the loss
        loss = criterion(out, edge_label)
        # backpropagate the gradients
        loss.backward()
        # update the model parameters
        optimizer.step()
        # evaluate the validation AUC
        val_auc = eval_link_predictor(model, val_data)
        test_auc = eval_link_predictor(model, test_data)
        # Store the valudation AUC and check if it has increased
        validation_auc_list.append(val_auc)
        last_10_auc_mean = sum(validation_auc_list[-10:]) / 10.0
        if epoch >10:
            validation_auc_list.pop(0)
        # Check if the validation auc has improved
        if last_10_auc_mean > best_val_auc:
            # update the best validation AUC score
            if epoch >= 1000:
                best_val_auc = last_10_auc_mean
            # reset the number of epochs without improvement
            num_epoche_no_improve = 0
            # save the best model as a pickle file
            with open(fname_output + '_model.pkl', 'wb') as f:
                pickle.dump(model, f)
            with open(fname_output + '.log', 'a') as f:
                f.write('..model saved\n')
                log_text = f"Epoch: {epoch:03d}, Train Loss: {loss:.3f}, Best Val AUC: {best_val_auc:.3f}, Val AUC: {val_auc:.3f}, Test AUC: {test_auc:.3f}\n"
                f.write(log_text)
        else:
            # increment the number of epochs without improvement
            num_epoche_no_improve += 1
        if epoch % 10 == 0:
            with open(fname_output + '.log', 'a') as f:
                log_text = f"Epoch: {epoch:03d}, Train Loss: {loss:.3f}, Best Val AUC: {best_val_auc:.3f}, Val AUC: {val_auc:.3f}, Test AUC: {test_auc:.3f}\n"
                f.write(log_text)
    return model


def get_network(path_files, TEST_ID, rng_seed, neg_ratio, file_suffix=''):
    fname_output = 'GAT/data/'+TEST_ID+'_rng'+str(rng_seed)+'_neg'+str(neg_ratio) + file_suffix
    graph_for_predicting = build_graph(path_files, TEST_ID, 'predicting')
    dt_all_possible_edge_index = pd.read_csv(path_files+'/'+TEST_ID+'_edge_index_for_predicting.txt',
                                             header=None, sep='\t')
    torch_all_possible_edge_index = torch.from_numpy(
        np.array(dt_all_possible_edge_index).T).cuda()
    with open(fname_output + '_model.pkl', 'rb') as f:
        best_model = pickle.load(f)
    best_model.eval()
    z = best_model.encode(graph_for_predicting.x,
                          graph_for_predicting.edge_index)
    out = best_model.decode(
        z, torch_all_possible_edge_index).view(-1).sigmoid()
    prediction_scores = out.cpu().detach()
    dt_all_possible_edge_labels = pd.read_csv(path_files+'/'+TEST_ID+'_all_possible_edges.txt',
                                              header=None, sep='\t')
    dt_all_possible_edge_labels['scores'] = prediction_scores
    dt_all_possible_edge_labels.to_csv(
        fname_output + '_prediction_score.txt', sep='\t', header=None, index=False)

def plot_histogram(tensor_data):
    # plot_histogram(prediction_scores)
    data = tensor_data.numpy()  # Convert tensor to NumPy array
    plt.hist(data, bins='auto', alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.5)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Histogram')
    plt.show()

def convert_to_networkx(graph, n_sample=None):
    g = to_networkx(graph, node_attrs=["x"])
    y = graph.y.numpy()
    if n_sample is not None:
        sampled_nodes = random.sample(g.nodes, n_sample)
        g = g.subgraph(sampled_nodes)
        y = y[sampled_nodes]
    return g, y


def plot_graph(g, y):
    plt.figure(figsize=(9, 7))
    networkx.draw_spring(g, node_size=30, arrows=False, node_color=y)
    plt.show()


def eval_link_predictor(model, data):
    model.eval()
    z = model.encode(data.x, data.edge_index)
    out = model.decode(z, data.edge_label_index).view(-1).sigmoid()
    return roc_auc_score(data.edge_label.cpu().numpy(), out.cpu().detach())
