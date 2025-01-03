# %%
# imports up top, babies
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import pandas as pd
import torch
import torchvision
from tqdm import tqdm

# %%
# Argle parser
parser = ArgumentParser(
    prog='calc_alexnet_activations',
    description='Get PyTorch AlexNet predictions for a directory of videos.',
    formatter_class=ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    '-i',
    '--in_paths',
    nargs='*',
    type=str,
    help="Paths to videos. Input multiple paths for multiple videos"
)
parser.add_argument(
    '-o',
    '--out_path',
    type=str,
    help="Path to write out overall activations file for all videos input"
)

args = vars(parser.parse_args())

video_paths = args['in_paths']
out_path = args['out_path']

# %%
# Fire up AlexNet from torchvision pretrained
model_torch = torchvision.models.alexnet(weights='IMAGENET1K_V1')
preprocess = torchvision.models.AlexNet_Weights.IMAGENET1K_V1.transforms()

# Turn off backward gradients for all runs of the model
# Right now we just be inferencing
for param in model_torch.parameters():
    param.requires_grad = False
# I think doing this is supposed to also turn the backward gradients off.
# But I'm not confident that it does, so, this is hocus pocus to do both but can't hurt
model_torch.eval()

# %%
# Workhorse classifier doer
# Loops directly over a set of paths
# Doesn't use a Dataset object bc we need input file flexibility

preds_all = []
    
for in_file in tqdm(video_paths, desc='Stimulus videos'):

    # read in video as stack of image tensors
    video = torchvision.io.read_video(
        in_file,
        pts_unit='sec'
    )
    # None of the videos have audio, so discard that from the loaded tuple
    # Also for convenience, discard dict labeling fps so that the videos look like 4D imgs
    # with dims frames x channels x height x width ... which is NOT the default order!
    frames = video[0].permute((0, 3, 1, 2))

    frames = preprocess(frames)

    # Actually calculate stuff
    # as written, this only gets the final classifier probabilities
    preds = model_torch(frames)          

    # format predictions as dataframe for writing to csv at the end
    preds = pd.DataFrame(
        preds.numpy(),
        index = range(preds.size()[0])
    )
    preds.index.name = 'frame'
    # Manually add the video name on as an index now as well
    preds['video'] = os.path.split(in_file)[1]
    preds = preds.set_index('video', append=True)

    # Now bind onto the superlists
    preds_all.append(preds)

preds_all = pd.concat(preds_all)
preds_all.to_csv(out_path)
