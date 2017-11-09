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

imgB_R = imgB(:,:,1);
imgB_B = imgB(:,:,2);
imgB_G = imgB(:,:,3);

N = size(imgB); 
sizeX = N(2);
sizeY = N(1);
N = N(1) * N(2);

% Initialize cluster centers (equally distribute on the image). Each cluster is represented by 5D feature (L, a, b, x, y)
% Hint: use linspace, meshgrid
S = sqrt(N/K);
C = zeros(5,K,'double');     % setup as [x1 ... xn; y1 ... yn; r; b; g;]
Cavg = zeros(5,K,'double');  % setup as [x1 ... xn; y1 ... yn; r; b; g;]

last = [S; S];      % last = [x,y]
C(1:2,1) = last;
C(3,1) = interp2(imgB_R, last(1), last(2), 'bilinear');
C(4,1) = interp2(imgB_B, last(1), last(2), 'bilinear');
C(5,1) = interp2(imgB_G, last(1), last(2), 'bilinear');

for i = 2:K
    last(1) = last(1) + S;
    if last(1) >= sizeX
        last(1) = S/2;
        last(2) = last(2) + S;
    end
    C(1:2,i) = last;

    % compute the color seperately for each channel
    C(3,i) = interp2(imgB_R, last(1), last(2), 'bilinear');
    C(4,i) = interp2(imgB_B, last(1), last(2), 'bilinear');
    C(5,i) = interp2(imgB_G, last(1), last(2), 'bilinear');
    
    if isnan(C(3,i))
        disp('uh oh');
    end
    
    
%     % display on imgB
%     x = round(C(1,i)); y = round(C(2,i));
%     imgB(y,x,:) = [0,1,0];
end


% assign each pixel to the nearest cluster
% cMap stored as follows: y,x indicates pixel position: and then:
% [clust,dist]
cMap = ones(sizeY,sizeX,2);
cMap = -1 * cMap;   % everything not initialized set to 1
cCount = zeros(K);
m = 0.1;      % m constant

xprange = [0 0];
yprange = [0 0];

for cluster = 1:K         % for every cluster region, try to assign pixel
    cx = C(1,cluster);
    cy = C(2,cluster);
    xprange = [ceil(max(1,cx-S)), floor(min(sizeX,cx+S))];
    yprange = [ceil(max(1,cy-S)), floor(min(sizeY,cy+S))];
    
    % initially set the pixels in the search range
    for x=xprange(1):xprange(2)
        for y=yprange(1):yprange(2)
            p = [x,y,imgB(y,x,1),imgB(y,x,2),imgB(y,x,3)]';
            
            dc = sqrt(sum((p(3:5) - C(3:5,cluster)).^2)); % colors only
            ds = sqrt(sum((p(1:2) - C(1:2,cluster)).^2)); % spatial distance only
            d = sqrt(dc^2 + (ds/S)^2 * m^2);

            if cMap(y,x,1) == -1    % if not set yet
                cMap(y,x,1) = cluster;
                cMap(y,x,2) = d;
                
            else                    % if already been set
                if cMap(y,x,2) > d  % if previous distance is larger, assign this one
                    cMap(y,x,1) = cluster;
                    cMap(y,x,2) = d;
                    
                end
            end
            
        end
    end
end

% attempt to visualize the initial cluster
% colors = rand(K,3);    % randomly assign color to each cluster: k clusters, 3 color channels
% for x=1:sizeX
%     for y = 1:sizeY-2
%         imgB(y,x,:) = colors(cMap(y,x,1),:);
%     end
% end





figure(); imshow(imgB);
% a = 1;





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