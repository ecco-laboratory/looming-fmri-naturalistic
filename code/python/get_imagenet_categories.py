# %% 
# imports
# this notebook script expects monica's torchtrain environment
import sys
sys.path.append("/home/mthieu/Repos/emonet-py/python/myutils")
from emonet_utils import Cowen2017Dataset

import torch
from torchvision.models import resnet50, ResNet50_Weights
these_weights = ResNet50_Weights.IMAGENET1K_V2
import pandas as pd
from tqdm import tqdm

video_path = "/home/mthieu/Repos/emonet-py/ignore/datasets/subjective/stimuli/fps10"
metadata_path = "/home/mthieu/Repos/emonet-py/ignore/datasets/subjective/metadata/kragel2019_all_video_10fps_ids.csv"
# %% 
# init model
resnet = resnet50(weights = these_weights)
resnet.eval()

# %%
# read those motherfuckin videos in
# ideally these would have been all in one dataset bc for our purposes we don't need to obey the train/test split
# but I can't be arsed to rewrite the dataset method
ck_torchdata_train = Cowen2017Dataset(
    root=video_path,
    annPath=metadata_path,
    censor=True,
    train=True,
    transform=these_weights.transforms()
)

ck_torchdata_test = Cowen2017Dataset(
    root=video_path,
    annPath=metadata_path,
    censor=True,
    train=False,
    transform=these_weights.transforms()
)

# Set batch_size here to 1 so it's just one video at a time
# BUT! Each video effectively acts as a batch of frames, as long as time is in the first dim
ck_torchloader_train = torch.utils.data.DataLoader(ck_torchdata_train, batch_size=1)
ck_torchloader_test = torch.utils.data.DataLoader(ck_torchdata_test, batch_size=1)

# %%
ids_all = []
emotions_all = []
categories_all = []

for loader in [ck_torchdata_train, ck_torchdata_test]:
    for vid, lab in tqdm(loader):
        vid = vid.squeeze()
        ids_all.append(lab['id'])
        emotions_all.append(lab['emotion'])
        pred = resnet(vid)
        # average preds across frames of each video
        pred = pred.mean(axis=0)
        class_id = pred.argmax().item()
        category_name = these_weights.meta["categories"][class_id]
        categories_all.append(category_name)

# %%
# convert the frames x emotions predictions into dataframes
# and concatenate for writing out

df = pd.DataFrame({'emotion': emotions_all, 'object': categories_all}, index=ids_all)

df.to_csv("ignore/data/norm/ck2017_imagenet_categories.csv")
