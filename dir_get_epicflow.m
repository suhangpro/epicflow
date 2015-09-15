function dir_get_epicflow(dir_path, save_path, n_pool, opt_step, opt_interval)
%DIR_GET_EPICFLOW Compute EpicFlow for frames in a folder 
% dir_path 
%   directory containing consecutive frames (.png files)
% save_path: (default: folder_path)
%   computed flow files will be saved here
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
root_path = fileparts(mfilename('fullpath'));
addpath(genpath(root_path));

% input
if exist(dir_path,'dir')~=7, error('Wrong path: %s', dir_path); end
if ~exist('save_path','var') || isempty(save_path), 
    save_path = dir_path;
end;
if ~exist('n_pool','var'), n_pool = 0; end
if ~exist('opt_step','var'), opt_step = 1; end
if ~exist('opt_interval','var'), opt_interval = 1; end

update_dir(dir_path, save_path, n_pool, opt_step, opt_interval); 

end

function update_dir( cur_dir, cur_save_dir, n_pool, opt_step, opt_interval)

% frame & flow files 
files = dir(fullfile(cur_dir,'*.png')); 
file_names = sort({files.name}); 
dirs = dir(cur_dir);
dir_names = {dirs.name};
dir_names = setdiff(dir_names(cell2mat({dirs.isdir})),{'.','..'});

if ~exist(cur_save_dir,'dir'), mkdir(cur_save_dir); end;

% do work for files
if ~isempty(file_names),
    fprintf('Computing optical flows for %s ...\n', cur_dir); 
    
    first_idxs = 1:opt_step:(numel(file_names)-opt_interval);
    second_idxs = first_idxs + opt_interval;
    first_files = cellfun(@(s) fullfile(cur_dir, s), file_names(first_idxs), ...
        'UniformOutput', false);
    second_files = cellfun(@(s) fullfile(cur_dir, s), file_names(second_idxs), ...
        'UniformOutput', false);
    flo_files = cellfun(@(s) fullfile(cur_save_dir, [s(1:end-4) '.flo']), ...
        file_names(first_idxs), 'UniformOutput', false);
    
    % get epicflow
    N = numel(flo_files);
    if n_pool>0,
        pool = gcp('nocreate');
        if isempty(pool),
            parpool(n_pool);
        elseif pool.NumWorkers<n_pool,
            delete(pool);
            parpool(n_pool);
        end
        parfor_progress(N);
        parfor i=1:N,
            get_epicflow(first_files{i}, second_files{i}, flo_files{i}, []);
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
end

% process subdirs recursively
for d=1:numel(dir_names), 
    update_dir(fullfile(cur_dir,dir_names{d}), ...
        fullfile(cur_save_dir,dir_names{d}), n_pool, opt_step, opt_interval); 
end

end
