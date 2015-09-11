function f = get_epicflow( im1_path, im2_path, flow_path, model_path)
%GET_EPICFLOW Wrapper for epicflow 
%   EpicFlow: Edge-Preserving Interpolation of Correspondences for Optical 
%   Flow, Jerome Revaud, Philippe Weinzaepfel, Zaid Harchaoui and Cordelia 
%   Schmid, CVPR 2015.

% setup 
dir_path = fileparts(mfilename('fullpath'));
SED_path = fullfile(dir_path,'SED_edge');
piotr_path = fullfile(dir_path,'piotr_toolbox');
dm_path = fullfile(dir_path,'deepmatching_1.2_c++','deepmatching-static');
epicflow_path = fullfile(dir_path,'EpicFlow_v1.00','epicflow-static');

addpath(genpath(SED_path));
addpath(genpath(piotr_path));

processid = feature('getpid');

% input
if ischar(im1_path), 
    if im1_path(1)~='/', im1_path = fullfile(pwd,im1_path); end;
    im1 = imread(im1_path); 
else
    im1 = im1_path;
    im1_path = fullfile(pwd,sprintf('im1_%d.tmp.png',processid));
    imwrite(im1,im1_path);
end
if ischar(im2_path), 
    if im2_path(1)~='/', im2_path = fullfile(pwd,im2_path); end;
else
    im2 = im2_path;
    im2_path = fullfile(pwd,sprintf('im2_%d.tmp.png',processid));
    imwrite(im2,im2_path);
end
if size(im1,3)==1, im1 = repmat(im1,[1 1 3]); end;

if ~exist('model_path','var') || isempty(model_path), 
    model_path = fullfile(fileparts(epicflow_path),'modelFinal.mat');
end
if ischar(model_path), 
    load(model_path); 
    if ~exist('model','var'), error('model not found'); end
else
    model = model_path;
end

if flow_path(1)~='/', flow_path = fullfile(pwd,flow_path); end

% extract edges 
edges = edgesDetect(im1,model);
edge_path = fullfile(pwd,sprintf('edge_%d.tmp',processid));
fid = fopen(edge_path,'wb');
fwrite(fid,transpose(edges),'single');
fclose(fid);

% matching
match_path = fullfile(pwd,sprintf('match_%d.tmp',processid));
cmd = sprintf('%s %s %s -png_settings -out %s', dm_path, ...
    im1_path, im2_path, match_path);
if(system(cmd)), error('Error while executing: %s',cmd); end

% epicflow
cmd = sprintf('%s %s %s %s %s %s', epicflow_path, im1_path, im2_path, ...
    edge_path, match_path, flow_path); 
if(system(cmd)), error('Error while executing: %s',cmd); end
f = flow_path; 

% clean up
if ~isempty(strfind(im1_path,int2str(processid))), delete(im1_path); end
if ~isempty(strfind(im2_path,int2str(processid))), delete(im2_path); end
delete(edge_path);
delete(match_path);

end

