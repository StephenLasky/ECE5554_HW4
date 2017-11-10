% Add path
addpath(genpath('GCmex1.5'));
im = im2double( imread('cat.jpg') );

org_im = im;
im = rgb2gray(im);

H = size(im, 1); W = size(im, 2); K = 3;
num_components = 5;

% from piazza:
% You can use the boolean matrix we provided to initialize the label, 
% then you solve a K-component GMM model (fit one GMM for foreground region 
% and one GMM for background). After you get the model, you can calculate the 
% probability of being foreground/background of each pixel, using the simple 
% formula of Gaussian distribution.

% Load the mask
load cat_poly
inbox = poly2mask(poly(:,1), poly(:,2), size(im, 1), size(im,2));

% 1) Fit Gaussian mixture model for foreground regions
% mu_fg = mean(im(inbox));
% sigma_fg = sqrt(mean((im(inbox)-mu_fg).^2))
GMM_FG = fitgmdist(im(inbox),num_components);


% 2) Fit Gaussian mixture model for background regions
% mu_bg = mean(im(~inbox))
% sigma_bg= sqrt(mean((im(~inbox)-mu_bg).^2))
GMM_BG = fitgmdist(im(~inbox),num_components);

% 3) Prepare the data cost
% - data [Height x Width x 2] 
% - data(:,:,1) the cost of assigning pixels to label 1
% - data(:,:,2) the cost of assigning pixels to label 2
data = zeros(H,W,2,'double');
% cost_FG = - log(pdf(GMM_FG, X));
% cost_BG = - log(pdf(GMM_BG, X));

for i = 1:H
    data(i,:,1) = pdf(GMM_FG,im(i,:)')';
    data(i,:,2) = pdf(GMM_BG,im(i,:)')';
end

data = arrayfun(@(x) -log(x), data);
data_range = [min(min(min(data))), max(max(max(data)))]



% 4) Prepare smoothness cost
% - smoothcost [2 x 2]
% - smoothcost(1, 2) = smoothcost(2,1) => the cost if neighboring pixels do not have the same label
smoothcost = [0 1; 1 0];
% smoothcost = [0.1 0.5 ; 0.5 0.1];

% 5) Prepare contrast sensitive cost
% - vC: [Height x Width]: vC = 2-exp(-gy/(2*sigma)); 
% - hC: [Height x Width]: hC = 2-exp(-gx/(2*sigma));
sigma = 1;
[gx, gy] = imgradientxy(im);
vC = 2-exp(-gy/(2*sigma)); 
hC = 2-exp(-gx/(2*sigma));

vC_range = [min(min(vC)), max(max(vC))]
hC_range = [min(min(hC)), max(max(hC))]

% 6) Solve the labeling using graph cut
% - Check the function GraphCut
[gch] = GraphCut('open', data ,smoothcost, vC, hC);
[gch labels] = GraphCut('expand', gch);
size(gch)
size(labels)


% 7) Visualize the results
% im(:,:,1) = im(:,:,1) .* inbox;
% im(:,:,2) = im(:,:,2) .* inbox;   % for G
% im(:,:,3) = im(:,:,3) .* inbox;   % for B
% imshow(im);