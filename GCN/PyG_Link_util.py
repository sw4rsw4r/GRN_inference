from networkx import *
import torch
import torch_geometric.transforms as T
from torch_geometric.nn import GCNConv
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
random.seed(123)


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
    graph = Data(edge_index=edge_index)

    embeddings = load_embedding(path_files, TEST_ID, type)
    graph.num_nodes = len(embeddings)
    graph.x = torch.from_numpy(embeddings).type(torch.float32)
    return (graph)


class Net(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, out_channels)

    def encode(self, x, edge_index):
        x = self.conv1(x, edge_index).relu()
        return self.conv2(x, edge_index)

    def decode(self, z, edge_label_index):
        return (z[edge_label_index[0]] * z[edge_label_index[1]]).sum(
            dim=-1
        )  # product of a pair of nodes on each edge

    def decode_all(self, z):
        prob_adj = z @ z.t()
        return (prob_adj > 0).nonzero(as_tuple=False).t()


def train_link_predictor(
    model, train_data, val_data, test_data, optimizer, criterion, TEST_ID, rng_seed=1234, neg_ratio=1.0, file_suffix=''):
    fname_output = 'GCN/data/'+TEST_ID+'_rng'+str(rng_seed)+'_neg'+str(neg_ratio) + file_suffix
    best_val_auc = 0.
    validation_auc_list = []
    patience = 500.
    num_epoche_no_improve = 0
    with open( fname_output + '.log', 'w') as f:
        f.write('Start training..\n')
    n_epochs = 5000
    for epoch in range(1, n_epochs + 1):
        model.train()
        optimizer.zero_grad()
        z = model.encode(train_data.x, train_data.edge_index)
        # sampling training negatives for every training epoch
        neg_edge_index = negative_sampling(
            edge_index=train_data.edge_index, num_nodes=train_data.num_nodes,
            num_neg_samples=train_data.edge_label_index.size(1), method='sparse')
        edge_label_index = torch.cat(
            [train_data.edge_label_index, neg_edge_index],
            dim=-1,
        )
        edge_label = torch.cat([
            train_data.edge_label,
            train_data.edge_label.new_zeros(neg_edge_index.size(1))
        ], dim=0)
        out = model.decode(z, edge_label_index).view(-1)
        loss = criterion(out, edge_label)
        loss.backward()
        optimizer.step()
        val_auc = eval_link_predictor(model, val_data)
        test_auc = eval_link_predictor(model, test_data)
        validation_auc_list.append(val_auc)
        last_10_auc_mean = sum(validation_auc_list[-10:]) / 10.0
        if epoch >10:
            validation_auc_list.pop(0)
        if last_10_auc_mean > best_val_auc:
            if epoch >= 1000:
                best_val_auc = last_10_auc_mean
            num_epoche_no_improve = 0
            with open(fname_output + '_model.pkl', 'wb') as f:
                pickle.dump(model, f)
            with open(fname_output + '.log', 'a') as f:
                f.write('..model saved\n')
                log_text = f"Epoch: {epoch:03d}, Train Loss: {loss:.3f}, Best Val AUC: {best_val_auc:.3f}, Val AUC: {val_auc:.3f}, Test AUC: {test_auc:.3f}\n"
                f.write(log_text)
        else:
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
    z = best_model.encode(graph_for_predicting.x, graph_for_predicting.edge_index)
    out = best_model.decode(z, torch_all_possible_edge_index).view(-1).sigmoid()
    prediction_scores = out.cpu().detach()
    dt_all_possible_edge_labels = pd.read_csv(path_files+'/'+TEST_ID+'_all_possible_edges.txt',
                                              header=None, sep='\t')
    dt_all_possible_edge_labels['scores'] = prediction_scores
    dt_all_possible_edge_labels.to_csv(
            fname_output + '_prediction_score.txt', sep='\t', header=None, index=False)


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
