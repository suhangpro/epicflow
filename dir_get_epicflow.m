function dir_get_epicflow( folder_path, flow_path, model_path, n_pool, opt_step, opt_interval)
%DIR_GET_EPICFLOW Compute EpicFlow for frames in a folder 
% folder_path 
%   path to a folder containing consecutive frames (.png files)
% flow_path: (default: folder_path)
%   computed flow files will be saved here
% model_path: (default: '')
%   edge detection model
% n_pool: (default:: 0) 
%   number of cpu cores to use
% opt_step: (default:: 1)
%   the 1st frames will be OPT_STEP frames away
% opt_interval: (default: 1)
%   the two frames will be OPT_INTERVAL frames away
% 
%   EpicFlow: Edge-Preserving Interpolation of Correspondences for Optical 
%   Flow, Jerome Revaud, Philippe Weinzaepfel, Zaid Harchaoui and Cordelia 
%   Schmid, CVPR 2015.
% 

% setup 
dir_path = fileparts(mfilename('fullpath'));
addpath(genpath(dir_path));

% input
if exist(folder_path,'dir')~=7, error('Wrong path: %s', folder_path); end
if ~exist('flow_path','var') || isempty(flow_path), 
    flow_path = folder_path;
elseif ~exist(flow_path,'dir'), mkdir(flow_path); end;
if ~exist('model_path','var'), model_path = []; end 
if ~exist('n_pool','var'), n_pool = 0; end
if ~exist('opt_step','var'), opt_step = 1; end
if ~exist('opt_interval','var'), opt_interval = 1; end

% frame & flow files 
files = dir(fullfile(folder_path,'*.png')); 
files = {files.name}; 
first_idxs = 1:opt_step:(numel(files)-opt_interval); 
second_idxs = first_idxs + opt_interval; 
first_files = cellfun(@(s) fullfile(folder_path, s), files(first_idxs), ...
    'UniformOutput', false);
second_files = cellfun(@(s) fullfile(folder_path, s), files(second_idxs), ...
    'UniformOutput', false);
flo_files = cellfun(@(s) fullfile(flow_path, [s(1:end-4) '.flo']), ...
    files(first_idxs), 'UniformOutput', false);  

% get epicflow
N = numel(flo_files); 
if n_pool>0, 
    pool = parpool(n_pool); 
    parfor_progress(N);
    parfor i=1:N, 
        get_epicflow(first_files{i}, second_files{i}, flo_files{i}, model_path);
        parfor_progress();
    end
else
    parfor_progress(N); 
    for i=1:N,
        get_epicflow(first_files{i}, second_files{i}, flo_files{i}, model_path);
        parfor_progress();
    end
end
parfor_progress(0);


