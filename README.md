# epicflow

A Matlab wrapper for EpicFlow:

EpicFlow: Edge-Preserving Interpolation of Correspondences for Optical 
Flow, Jerome Revaud, Philippe Weinzaepfel, Zaid Harchaoui and Cordelia 
Schmid, CVPR 2015.

Tested only on 64 bit Linux. 

## Usage 

* compute optical flow
```
>> f = get_epicflow(IMAGE1_PATH, IMAGE2_PATH, FLOW_SAVE_PATH);
```
* compute optical flow for a folder (see dir_get_epicflow.m for more options)
```
>> dir_get_epicflow(FOLDER_PATH, FLOW_SAVE_FOLDER_PATH)
```
* visualize result (uses [Middlebury optical flow toolbox](http://vision.middlebury.edu/flow/code/flow-code-matlab.zip))
```
>> addpath(genpath('utils'));
>> imshow(flowToColor(readFlowFile(f)));
```
