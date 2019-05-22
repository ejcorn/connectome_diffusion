function cropped_image = crop_brain_image(cdata,bg,marg)

% cdata = X-by-Y-by-3 (RGB) image
% bg = value of background. Should be 255, but for some reason it is 
% 240 in the case of frame2im output
% marg = margin of background to leave as a fraction of original image size
% take in RGB image of two brain surface plots next to each other
% looks like [X X; O O] where X is brain and O is white space
% crop images to look like [X X]
if ~exist('marg','var')
    marg = 0.05;
end

if ~exist('bg','var')
    bg = 240;
end

X = sum(cdata,3); % average over colors
[x_sz,y_sz] = size(X);
[x,y] = find(X ~= 3*bg);
crop_top = min(x) - round(marg*x_sz);
crop_bottom = max(x) + round(marg*x_sz);
crop_left = min(y) - round(marg*y_sz);
crop_right = max(y) + round(marg*y_sz);
cropped_image = cdata(crop_top:crop_bottom,crop_left:crop_right,:);


