# Human superior colliculus pathways represent the form and motion of looming objects
The superior colliculus is well known for its roles in visual orienting, oculomotor control, attention, and defensive behavior across species. Recently, [we predicted and found](https://www.cell.com/iscience/fulltext/S2589-0042(24)01108-8) that representations from a shallow convolutional neural network could predict defensive blinking to looming objects in infants and superior colliculus responses to optical expansion in adults. These findings suggest that the superior colliculus may coordinate defensive responses to looming in humans. 

In this project, we tested whether the human superior colliculus functions in isolation during looming threat perception, or if it covaries with cortical networks involved in visual salience and object recognition. We used computational models of [looming detection](https://elifesciences.org/articles/72067), [visual saliency](https://ieeexplore.ieee.org/document/730558/), and [object recognition](https://papers.nips.cc/paper_files/paper/2012/hash/c399862d3b9d6b76c8436e924a68c45b-Abstract.html) to predict patterns of superior colliculus BOLD response acquired as participants viewed a series of naturalistic videos or performed a working memory task, and to examine their functional connectivity with the rest of the brain. A detailed description of this work is forthcoming. This repository provides the source code used to conduct these analyses.

<p align="center">
<img src="https://github.com/ecco-laboratory/looming-fmri-naturalistic/blob/master/images/GraphicalAbstract.png" width="800">
</p>


## Dependencies 
This code uses [Canlab Core Tools](https://github.com/canlab/CanlabCore/tree/master), which is an object oriented toolbox that uses the [SPM software](https://www.fil.ion.ucl.ac.uk/spm/) to process fMRI data. Instructions for installing CANLab Core Tools which provided many of the functions and Neuroimaging Pattern Masks used in the analyses can be found [here](https://canlab.github.io/_pages/canlab_help_1_installing_tools/canlab_help_1_installing_tools.html).

Instructions for installing SPM12 can be found [here](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/).

The code requires both [PyTorch](https://pytorch.org/) to run. Installation should take approximately 5 minutes.

Analyses were performed using MATLAB R2024a on Ubuntu 20.04.6 LTS (Focal Fossa). No non-standard or specialized hardware was required. 
