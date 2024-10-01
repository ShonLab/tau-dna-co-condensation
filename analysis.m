nframe = 151;
mov_DNA = loadTifStack16('raw\Single DNA dororok_DNA.tif',nframe);
mov_tau = loadTifStack16('raw\Single DNA dororok_tau.tif',nframe);

% registration
img_avg_DNA = rescale(mean(mov_DNA(:,:,end-10:end),3),0,5);
img_avg_tau = rescale(mean(mov_tau(:,:,end-10:end),3),0,5);

% [movPts,fxdPts] = cpselect(img_avg_tau,img_avg_DNA, 'Wait',true);
tform = fitgeotform2d(movPts,fxdPts,"affine");
img_avg_tau_reg = imwarp(img_avg_tau,tform,OutputView=imref2d(size(img_avg_tau)));

maxfig(1); clf;
subplot(211);
imshow2(cat(3,img_avg_tau,img_avg_DNA,zeros(size(img_avg_DNA))));
subplot(212);
imshow2(cat(3,img_avg_tau_reg,img_avg_DNA,zeros(size(img_avg_DNA))));

for f = 1:nframe
    mov_tau(:,:,f) = imwarp(mov_tau(:,:,f),tform,OutputView=imref2d(size(img_avg_tau)));
end
% save('analysis');

%%
hex_DNA = '#7FFF00';
hex_tau = '#CC79A7';
cmap_DNA = hex2map(hex_DNA);
cmap_tau = hex2map(hex_tau);

rrange = 188+(1:110); crange = 305+(1:150);
Imin_DNA = 120; Imax_DNA = 250;
Imin_tau = 106; Imax_tau = 118;

%% export movie
video = VideoWriter('Single DNA dororok_cropped.avi', 'Uncompressed AVI');
video.FrameRate = 10;
open(video);
for f = 1:nframe
    img_DNA = rescale(mov_DNA(rrange,crange,f), 0,255, 'InputMin',Imin_DNA,'InputMax',Imax_DNA);
    img_DNA_rgb = ind2rgb(round(img_DNA),cmap_DNA);
    img_tau = rescale(mov_tau(rrange,crange,f), 0,255, 'InputMin',Imin_tau,'InputMax',Imax_tau);
    img_tau_rgb = ind2rgb(round(img_tau),cmap_tau);
    img_merge = img_DNA_rgb + img_tau_rgb;
    img_tmp = [img_DNA_rgb, img_tau_rgb, img_merge];
    writeVideo(video, uint8(img_tmp*255));
end
close(video);

%%
maxfig(1); clf;
f_list = [21,41,61,81,101];
img_sel = zeros(numel(rrange),numel(crange),3,2,numel(f_list));
for fi = 1:5
    f = f_list(fi);
    img_DNA = rescale(mov_DNA(rrange,crange,f), 0,255, 'InputMin',Imin_DNA,'InputMax',Imax_DNA);
    img_DNA_rgb = ind2rgb(round(img_DNA),cmap_DNA);
    img_tau = rescale(mov_tau(rrange,crange,f), 0,255, 'InputMin',Imin_tau,'InputMax',Imax_tau);
    img_tau_rgb = ind2rgb(round(img_tau),cmap_tau);
    img_merge = img_DNA_rgb + img_tau_rgb;
    
    subplot(3,5,fi);
    imshow2(img_DNA_rgb);
    
    subplot(3,5,fi+5);
    imshow2(img_tau_rgb);

    subplot(3,5,fi+10);
    imshow2(img_merge);

    img_sel(:,:,:,1,fi) = img_DNA_rgb;
    img_sel(:,:,:,2,fi) = img_tau_rgb;
end

%% sample kymograph
figure(2); clf;
rrange_sub = 89+(-5:5); crange_sub = 77+(-10:65);

kym_DNA = squeeze(mean(mov_DNA(rrange(rrange_sub),crange(crange_sub),:),1));
kym_DNA = rescale(kym_DNA, 0,255, 'InputMin',Imin_DNA,'InputMax',Imax_DNA)';
kym_DNA_rgb = ind2rgb(round(kym_DNA),cmap_DNA);

kym_tau = squeeze(mean(mov_tau(rrange(rrange_sub),crange(crange_sub),:),1));
kym_tau = rescale(kym_tau, 0,255, 'InputMin',Imin_tau,'InputMax',Imax_tau)';
kym_tau_rgb = ind2rgb(round(kym_tau),cmap_tau);

imshow2(kym_DNA_rgb+kym_tau_rgb);

%% length change
frange = 122:151;
[nmol,xpos,ypos] = countSM(mean(mov_DNA(:,:,frange),3),100,0);
sel_n = setdiff(1:nmol,[1,3,4,6,8,10,18,19,21:23,25,27:32,34,36,37,38,40,44]);
xpos = xpos(sel_n);
ypos = ypos(sel_n);
nmol = numel(sel_n);

figure(3); clf;
img_DNA_bkg = double(min(mov_DNA,[],3));
img_tau_bkg = double(min(mov_tau,[],3));

subplot(211);
frange1 = 1:30;
imshow2(mean(mov_DNA(:,:,frange1),3)-img_DNA_bkg,[0,50]); hold all;
plot(xpos,ypos,'y.');
text(xpos,ypos,arrayfun(@num2str,1:nmol,'unif',0),'hori','right','color','y');

subplot(212);
imshow2(mean(mov_DNA(:,:,frange),3)-img_DNA_bkg,[0,50]);

maxfig(4); clf;
maxfig(5); clf;
length = zeros(nmol,2);
I_end = zeros(21,nmol,2);
for n = 1:nmol
    % r0 = 277; c0 = 382;
    r0 = round(ypos(n)); c0 = round(xpos(n));
    mov_tmp = double(mov_DNA(r0+(-5:5),c0+(-10:65),:));
    mov_tmp = mov_tmp - repmat(img_DNA_bkg(r0+(-5:5),c0+(-10:65)),[1,1,nframe]);
    kym_tmp = squeeze(mean(mov_tmp,1));

    I_tmp = max(kym_tmp(:,1:51),[],2);
    I_min = mean(kym_tmp(end-4:end,end-30:end),'all');
    I_max = mean(I_tmp(11:40));
    I_threshold = I_min+(I_max-I_min)*.7;
    if ismember(n,[11,14,15])
        yend = find(smooth(I_tmp(1:65),3)>I_threshold,1,'last');
    else
        yend = find(smooth(I_tmp,3)>I_threshold,1,'last');
    end
    length(n,1) = (yend-6)*.1083;

    I_tmp2 = mean(kym_tmp(:,end-10+(-1:1)),2);
    I_end(:,n,1) = I_tmp2(1:21);
    
    figure(4);
    subplot(3,7,n);
    imshow2(kym_tmp,[0,20]);
    % vline(31);
    hline(yend);
    title(n);

    figure(5);
    subplot(6,7,n);
    plot(I_tmp);
    ylim([0,200]);
    hline(I_min,'k');
    hline(I_threshold,'r');
    vline(yend,'b');

    subplot(6,7,n+21);
    plot(I_tmp2);
    y1 = interp1(I_tmp2(1:11),1:11,max((I_tmp2)/2));
    y2 = interp1(I_tmp2(11:21),11:21,max((I_tmp2)/2));
    % length(n,2) = max((y2-y1)/2*.1083-.532/2,0);
    length(n,2) = (y2-y1)/2*.1083-.532/2;

    % ylim([0,200]);
    % hline(I_min,'k');
    % hline(I_threshold,'r');
    % vline(yend,'b');

    mov_tmp = double(mov_tau(r0+(-5:5),c0+(-10:65),:));
    mov_tmp = mov_tmp - repmat(img_tau_bkg(r0+(-5:5),c0+(-10:65)),[1,1,nframe]);
    kym_tmp = squeeze(mean(mov_tmp,1));
    I_tmp2 = mean(kym_tmp(:,end-10+(-1:1)),2);
    I_end(:,n,2) = I_tmp2(1:21);    
end

figure(6); clf;
bar(mean(length),'facecolor',hex2rgb(hex_DNA)); hold all;
errorbar(mean(length),std(length),'k','linestyle','none');
plotSpread(length);

save('analysis');