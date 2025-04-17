from GAT_hetero.GAT_hetero_util import *
import os
os.chdir('/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/')
TEST_ID = 'TEST3'

# 1. Build graphs for training and predicting
path_files = '/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/input_data_processing/data'
graph_for_training = build_graph(path_files, TEST_ID, 'training')

# 2. Split the graph into training and validation
split = T.RandomLinkSplit(
    num_val=0.05,
    num_test=0.,
    is_undirected=False,
    add_negative_train_samples=True,
    neg_sampling_ratio=1.0,
    edge_types=("TF", "interact", "gene")
)
train_data, val_data, test_data = split(graph_for_training)

# 3. Build a model and train
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = Net(train_data.x.shape[1], 128, 64).to(device)

model = to_hetero(model, train_data.metadata(), aggr = 'sum')

optimizer = torch.optim.Adam(params=model.parameters(), lr=0.002)
criterion = torch.nn.BCEWithLogitsLoss()
model = train_link_predictor(
    model, train_data, val_data, optimizer, criterion, n_epochs=5000)

# 4. Save the final predicted model to a file
get_network(path_files, TEST_ID, model)
