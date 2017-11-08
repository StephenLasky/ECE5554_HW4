% Add path
addpath(genpath('GCmex1.5'));
im = im2double( imread('cat.jpg') );

org_im = im;

H = size(im, 1); W = size(im, 2); K = 3;

% Load the mask
load cat_poly
inbox = poly2mask(poly(:,1), poly(:,2), size(im, 1), size(im,2));

% 1) Fit Gaussian mixture model for foreground regions

% 2) Fit Gaussian mixture model for background regions

% 3) Prepare the data cost
% - data [Height x Width x 2] 
% - data(:,:,1) the cost of assigning pixels to label 1
% - data(:,:,2) the cost of assigning pixels to label 2

% 4) Prepare smoothness cost
% - smoothcost [2 x 2]
% - smoothcost(1, 2) = smoothcost(2,1) => the cost if neighboring pixels do not have the same label

% 5) Prepare contrast sensitive cost
% - vC: [Height x Width]: vC = 2-exp(-gy/(2*sigma)); 
% - hC: [Height x Width]: hC = 2-exp(-gx/(2*sigma));

% 6) Solve the labeling using graph cut
% - Check the function GraphCut

% 7) Visualize the results