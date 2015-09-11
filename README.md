# epicflow

A Matlab wrapper for EpicFlow:

EpicFlow: Edge-Preserving Interpolation of Correspondences for Optical 
Flow, Jerome Revaud, Philippe Weinzaepfel, Zaid Harchaoui and Cordelia 
Schmid, CVPR 2015.

Tested only on 64 bit Linux. 

## Usage 

* compute optical flow
```
>> addpath(PATH_TO_REPO_ROOT);
>> f = get_epicflow(IMAGE1_PATH, IMAGE2_PATH, FLOW_SAVE_PATH);
```
* visualize result (uses [Middlebury optical flow toolbox](http://vision.middlebury.edu/flow/code/flow-code-matlab.zip))
```
>> imshow(flowToColor(readFlowFile(f)))
```
