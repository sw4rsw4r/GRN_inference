from GCN.PyG_Link_util import *
import os
os.chdir('/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/')
#TEST_ID = 'TEST3'
TEST_ID = 'NodeAb13_control_Mesoderm'
neg_ratio = 3.0


# 1. Build graphs for training and predicting
#path_files = '/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/input_data_processing/data'
path_files = '/home/seongwonhwang/Desktop/projects/git/Node_ablation_practice/GRN/data'
graph_for_training = build_graph(path_files, TEST_ID, 'training')

for rng_seed in range(1, 50 + 1):
    # Set a fixed seed for reproducibility
    random.seed(rng_seed)
    torch.manual_seed(rng_seed)
    np.random.seed(rng_seed)
    # 2. Split the graph into training and validation
    split = T.RandomLinkSplit(
        num_val=0.05,
        num_test=0.05,
        is_undirected=False,
        add_negative_train_samples=True,
        neg_sampling_ratio=neg_ratio,
    )
    train_data, val_data, test_data = split(graph_for_training)
    # 3. Build a model and train
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = Net(train_data.x.shape[1], 128, 64).to(device)
    optimizer = torch.optim.Adam(params=model.parameters(), lr=0.002)
    criterion = torch.nn.BCEWithLogitsLoss()
    model = train_link_predictor(
        model, train_data, val_data, test_data, optimizer, criterion, TEST_ID, rng_seed, neg_ratio)
    # 4. Save the final predicted model to a file
    #get_network(path_files, TEST_ID, rng_seed, neg_ratio)
