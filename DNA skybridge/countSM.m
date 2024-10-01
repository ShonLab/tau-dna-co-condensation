function [nmol,xpos,ypos] = countSM(img_avg,threshold,dispres)
% threshold = 7: good for single molecules. 

img_avg = double(img_avg);

% subtract background
% img_avg_filtered = medfilt2(img_avg,[21,21]);
img_avg_filtered = medfilt2(img_avg,[3,3]);
img_avg_bkgsub = img_avg-img_avg_filtered;

% filter image
H = fspecial('average',[3,3]);
img_avg_bkgsub = imfilter(img_avg_bkgsub,H);

% remove low-intensity peaks
img_avg_bkgsub(img_avg_bkgsub<threshold) = threshold;

% find molecules in the restricted area
Imol = imregionalmax(img_avg_bkgsub,8);
if all(all(Imol))
    Imol = false(size(img_avg_bkgsub));
end

Imask = true(size(img_avg));
Imask(6:end-5,6:end-5) = false;
Imol(Imask) = false;
[ypos,xpos] = find(Imol);
nmol = numel(ypos);

% verification by centroid
sel = false(nmol,1);
for k = 1:nmol
    subimg = img_avg(ypos(k)+(-2:2),xpos(k)+(-2:2));
    cnt = centroid(subimg);
    if sqrt(sum((cnt-3).^2)) < .7
        sel(k) = true;
        xpos(k) = xpos(k)+cnt(2)-3;
        ypos(k) = ypos(k)+cnt(1)-3;
    end
end
nmol = numel(find(sel));
xpos = xpos(sel);
ypos = ypos(sel);

% plot results
if dispres
    clf;
    ax1 = subplot(121);
    imshow3(img_avg,10);
    ax2 = subplot(122);
    imshow3(img_avg_bkgsub,10); hold on;
    plot(xpos,ypos,'yo','markersize',10);
    title(['nmol = ',num2str(nmol)]);
    linkaxes([ax1,ax2],'xy');
end