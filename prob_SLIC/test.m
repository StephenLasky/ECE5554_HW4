% import image
% im = im2double(imread('BSR/BSDS500/data/images/test/2018.jpg'));
% figure(); imshow(im);



% a = [1, 2, 3; 4, 5, 6; 7, 8, 9;];
% x = [1,1;2,2;3,3];
% diag(a(x(:,1),x(:,2)))


% % WORKS FOR 1-DIMENSIONAL COLOR
% color = [10, 20, 30];
% cMap = [1 2; 3 3];
% x = color(cMap)


% color = [10,15; 20,25; 30,35;];
% cMap = [1 2; 3 3];
% img = cMap;
% for r=1:2
%     img(r,:) = 
% end



cmap = ones(3,3,2);
cmap(:,:,1) = [1 2 3; 4 5 6; 7 8 9];
cmap(:,:,2) = [1 1 1; 2 1 1; 4 1 5];

d = [1 2 3; 0 5 6; 0 8 0];


S = d < cmap(:,:,1)

cmap(S) = d(S)





% A(A(:,:,2)>1) = 500

