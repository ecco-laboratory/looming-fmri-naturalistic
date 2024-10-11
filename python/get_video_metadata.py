# %%
# %%
# imports

import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from fractions import Fraction

import ffmpeg
import pandas as pd
from tqdm import tqdm

# %%
# Argle parser
parser = ArgumentParser(
    prog='get_video_metadata',
    description='Use ffprobe to extract metadata for a set of videos.',
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
    help="Path to write out overall metadata file for all videos input"
)

args = vars(parser.parse_args())

video_paths = args['in_paths']
out_path = args['out_path']


# %%
# Extract the metadata in a loop
out_all = {
    'width': [],
    'height': [],
    'duration': [],
    'frame_rate': [], 
    'frame_rate_float': [],
    'nframes': [],
    'video': []
}

for in_file in tqdm(video_paths, desc='Stimulus videos'):
    # USUALLY stream 0 is the video stream
    # but OCCASIONALLY? stream 1 is the video stream
    # so we need to check it
    metadata = ffmpeg.probe(in_file)['streams']
    if metadata[0]['codec_type'] == 'video':
        metadata = metadata[0]
    else:
        metadata = metadata[1]

    out_all['width'].append(metadata['width'])
    out_all['height'].append(metadata['height'])
    out_all['duration'].append(float(metadata['duration']))
    out_all['frame_rate'].append(metadata['r_frame_rate'])
    out_all['frame_rate_float'].append(float(Fraction(metadata['r_frame_rate'])))
    out_all['nframes'].append(float(metadata['nb_frames']))
    out_all['video'].append(os.path.split(in_file)[1])

out_all = pd.DataFrame(out_all)
out_all = out_all.set_index('video')
out_all.to_csv(out_path)
