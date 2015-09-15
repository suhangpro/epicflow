function opts = dir_resize( dir_path, save_path, opts )
% Resize images in the whole directory
% Hang Su
% 
% dir_path:     directory containing images, will be searched recursively 
% save_path:    (default:: dir_path) path to save resized images
% opts:    
%   .size       (default:: [256 256]) resizing target (NROWS,NCOLS)
%   .mode       ['crop']|'stretch'|'pad' 
%   .color      ['rgb']|'gray'
%   .ext        (default:: 'png') image extension
%   .pad_white  (default:: false) pad white instead of black


defaultOpts.size = [256 256]; % [height width]
defaultOpts.mode = 'crop'; 
defaultOpts.color = 'rgb'; 
defaultOpts.ext = 'png'; 
defaultOpts.pad_white = false; 

% return default options when called w/o any arguments
if nargin==0,
    opts = defaultOpts;
    return;
end

% if save_path is not specified, overwrite original image 
if ~exist('save_path','var') || isempty(save_path), 
  save_path = dir_path;
end

% populate the options not specified w/ default values
if ~exist('opts','var') || ~isstruct(opts), 
    opts = defaultOpts;
else
    missingFields = setdiff(fieldnames(defaultOpts),fieldnames(opts));
    for i=1:length(missingFields),
        opts.(missingFields{i}) = defaultOpts.(missingFields{i});
    end
end

update_dir(dir_path, save_path, opts); 

end

function update_dir(cur_dir, cur_save_dir, opts)

files = dir(fullfile(cur_dir, ['*.' opts.ext]));
file_names = {files.name}; 
dirs = dir(cur_dir);
dir_names = {dirs.name};
dir_names = setdiff(dir_names(cell2mat({dirs.isdir})),{'.','..'});

if~exist(cur_save_dir,'dir'), mkdir(cur_save_dir); end;

% do work
for i=1:numel(file_names), 
  im = imread(fullfile(cur_dir, file_names{i})); 
  im = im_resize(im,opts); 
  
  % save (possibly overwriting original image) 
  imwrite(im,fullfile(cur_save_dir, file_names{i})); 
  
end

for d = 1:numel(dir_names), 
  update_dir(fullfile(cur_dir,dir_names{d}), ...
    fullfile(cur_save_dir,dir_names{d}), opts);
end

end

function im = im_resize(im, opts)

if strcmpi(opts.color,'rgb'),   nch = 3;
elseif strcmpi(opts.color,'gray'), nch = 1;
else error('Color mode unsupported: %s',opts.color);
end
ar = opts.size(2)/opts.size(1); % aspect ratio

% convert color mode
if nch == 3 && size(im,3)==1,
    im = repmat(im,[1 1 3]);
elseif nch==1 && size(im,3)==3,
    im = rgb2gray(im);
end

% resize
sz = [size(im,1) size(im,2)];
if strcmpi(opts.mode,'crop'),         % mode 1: crop
    if sz(2)/sz(1) > ar,
        im_tmp = imresize(im,[opts.size(1) nan]);
        im = im_tmp(:,floor((size(im_tmp,2)-opts.size(2))/2)+(1:opts.size(2)),:);
    else
        im_tmp = imresize(im,[nan opts.size(2)]);
        im = im_tmp(floor((size(im_tmp,1)-opts.size(1))/2)+(1:opts.size(1)),:,:);
    end
elseif strcmpi(opts.mode,'stretch'),  % mode 2: stretch
    im = imresize(im, opts.size);
elseif strcmpi(opts.mode,'pad'),      % mode 3: pad
    if opts.pad_white,
        im_tmp = ones(opts.size(1),opts.size(2),nch,'uint8')*255;
    else
        im_tmp = zeros(opts.size(1),opts.size(2),nch,'uint8');
    end
    if sz(2)/sz(1) > ar,
        im = imresize(im,[nan opts.size(2)]);
        im_tmp(floor((opts.size(1)-size(im,1))/2)+(1:size(im,1)),:,:) = im;
    else
        im = imresize(im,[opts.size(1) nan]);
        im_tmp(:,floor((opts.size(2)-size(im,2))/2)+(1:size(im,2)),:) = im;
    end
    im = im_tmp;
end

end