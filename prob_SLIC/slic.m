function [cIndMap, time, imgVis] = slic(img, K, compactness)

%% Implementation of Simple Linear Iterative Clustering (SLIC)
%
% Input:
%   - img: input color image
%   - K:   number of clusters
%   - compactness: the weights for compactness
% Output: 
%   - cIndMap: a map of type uint16 storing the cluster memberships
%     cIndMap(i, j) = k => Pixel (i,j) belongs to k-th cluster
%   - time:    the time required for the computation
%   - imgVis:  the input image overlaid with the segmentation

% Put your SLIC implementation here

tic;
% Input data
imgB   = im2double(img);
cform  = makecform('srgb2lab');
imgLab = applycform(imgB, cform);

% Initialize cluster centers (equally distribute on the image). Each cluster is represented by 5D feature (L, a, b, x, y)
% Hint: use linspace, meshgrid


% SLIC superpixel segmentation
% In each iteration, we update the cluster assignment and then update the cluster center

numIter  = 10; % Number of iteration for running SLIC
for iter = 1: numIter
	% 1) Update the pixel to cluster assignment
	

	% 2) Update the cluster center by computing the mean

end

time = toc;


% Visualize mean color image
[gx, gy] = gradient(cIndMap);
bMap = (gx.^2 + gy.^2) > 0;
imgVis = img;
imgVis(cat(3, bMap, bMap, bMap)) = 1;

cIndMap = uint16(cIndMap);

end