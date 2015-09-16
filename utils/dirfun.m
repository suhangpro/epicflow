function dirfun( dir_path, save_path, func, file_pattern )
% dirfun Apply a function to each file in the directory
% 
% dir_path:     directory containing images, will be searched recursively 
% save_path:    (default:: dir_path) path to save resized images
% func          function that will be applied to each image found
% file_pattern: (default:: '*.png') images that will be processed
% 
% Hang Su

% if save_path is not specified, overwrite original image 
if ~exist('save_path','var') || isempty(save_path), 
  save_path = dir_path;
end

% default file pattern
if ~exist('file_pattern','var') || isempty(file_pattern), 
    file_pattern = '*';
end

% run recursively
update_dir(dir_path, save_path, func, file_pattern); 

end

function update_dir(cur_dir, cur_save_dir, func, file_pattern)

files = dir(fullfile(cur_dir, file_pattern));
file_names = {files.name}; 
dirs = dir(cur_dir);
dir_names = {dirs.name};
dir_names = setdiff(dir_names(cell2mat({dirs.isdir})),{'.','..'});

if~exist(cur_save_dir,'dir'), mkdir(cur_save_dir); end;

% do work
for i=1:numel(file_names), 
  im = imread(fullfile(cur_dir, file_names{i})); 
  im = func(im); 
  
  % save (possibly overwriting original image) 
  imwrite(im,fullfile(cur_save_dir, file_names{i})); 
  
end

for d = 1:numel(dir_names), 
  update_dir(fullfile(cur_dir,dir_names{d}), ...
    fullfile(cur_save_dir,dir_names{d}), func, file_pattern);
end

end
