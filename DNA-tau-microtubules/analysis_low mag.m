% %% load data
% % pth = 'raw3_low mag\';
% finfo = {dir([pth,'\*no tau_0*.tif']), dir([pth,'\*with tau_0*.tif'])};
% finfo{2}(end) = []; % ch 2 NO Tau + DNA + tubulin with tau_016.tif excluded
% npxr = 512; npxc = npxr;
% img = cell(2,1); nfile = zeros(2,1);
% for exp = 1:2
%     nfile(exp) = numel(finfo{exp});
%     img{exp} = zeros(npxr,npxc,nfile(exp),2);
%     for f = 1:nfile(exp)
%         disp([exp,f]);
%         % img_tmp = imread([pth,finfo{exp}(f).name]);
%         mov = loadTifStack16([pth,finfo{exp}(f).name]);
%         img{exp}(:,:,f,1) = mean(mov(1:512,256+(1:512),:),3);
%         img{exp}(:,:,f,2) = mean(mov(513:end,256+(1:512),:),3);
%     end
% end
% clear mov;
% 
% %% registration
% [movPts_all,fixPts_all] = deal([]);
% for f = 1:9
%     [movPts,fixPts] = cpselect(rescale(img{2}(:,:,f,2))*2,rescale(img{2}(:,:,f,1))*2,'Wait',true);
%     movPts_all = [movPts_all; movPts];
%     fixPts_all = [fixPts_all; fixPts];
% end
% tform = fitgeotform2d(movPts_all,fixPts_all,"affine");
% save('registration','*Pts_all','tform');

load('analysis_raw3');
load('registration');

%% before crosstalk correction
maxfig(11); clf;
Imax_list = [250,50];
% f_list = [5,5];
rrange = 412+(-99:100); crange = 226+(-99:100);
% rrange = 1:512; crange = 1:512;
img_bkg = medfilt2(imwarp(min(img{2}(:,:,:,2),[],3),tform,OutputView=imref2d([npxr,npxc])),[50,50]);
img_proc = zeros(numel(rrange),numel(crange),2,9,2);
for exp = 1:2
    for ch = 1:2
        for f = 1:9
            subplot(4,9,18*(exp-1)+9*(ch-1)+f);
            img_tmp = img{exp}(:,:,f,ch);
            if ch == 2
                img_tmp = imwarp(img_tmp,tform,OutputView=imref2d([npxr,npxc]));
                % img_tmp = img_tmp - img_bkg;
            end
            img_proc(:,:,exp,f,ch) = img_tmp(rrange,crange);
            img_tmp2 = rescale(img_tmp(rrange,crange),'InputMax',Imax_list(ch),'InputMin',0);
            imshow2(img_tmp2);
        end
    end
end

figure(12); clf;
tmp1 = img_proc(51:end-50,:,2,:,1);
tmp2 = img_proc(51:end-50,:,2,:,2);
sel = tmp2>250;
plot(tmp2(:),tmp1(:),'.'); hold all;
coeff = robustfit(tmp2(sel),tmp1(sel));
xdat = [min(tmp2(sel)),max(tmp2(sel))];
ydat = coeff(1) + coeff(2)*xdat;
plot(xdat,ydat,'r');

%% after crosstalk correction
maxfig(13); clf;
for exp = 1:2
    for ch = 1:2
        for f = 1:9
            subplot(4,9,18*(exp-1)+9*(ch-1)+f);
            img_tmp = img{exp}(:,:,f,ch);
            if ch == 1
                % img_cross = -100 + .6*imwarp(img{exp}(:,:,f,2),tform,OutputView=imref2d([npxr,npxc]));
                img_cross = coeff(1) + coeff(2)*imwarp(img{exp}(:,:,f,2),tform,OutputView=imref2d([npxr,npxc]));
                img_tmp = img_tmp - img_cross;
            elseif ch == 2
                img_tmp = imwarp(img_tmp,tform,OutputView=imref2d([npxr,npxc]));
                img_tmp = img_tmp - img_bkg;                
            end
            img_proc(:,:,exp,f,ch) = img_tmp(rrange,crange);
            img_tmp2 = rescale(img_tmp(rrange,crange));
            imshow2(img_tmp2);
        end
    end
end

%%
save('analysis_raw3');
save('analysis2','img_proc');

%% line detection
% title_list = {'without tau','with tau'};
% for exp = 1:2
%     for f = 1:9
%         img_tmp = imwarp(img{exp}(:,:,f,2),tform,OutputView=imref2d([npxr,npxc]));
%         img_tmp = img_tmp - img_bkg;
%         imwrite(uint16(img_tmp(101:end-50,51:end-50)),['tub img_',title_list{exp},'_',num2str(f),'.tif'])
%     end
% end

nline = zeros(2,9);
maxfig(1); clf;
for exp = 1:2
    maxfig(100+exp); clf;
    for f = 1:9
        figure(1); clf;
        img_tmp = imwarp(img{exp}(:,:,f,2),tform,OutputView=imref2d([npxr,npxc]));
        img_tmp = img_tmp - medfilt2(img_tmp,[20,20]);
        img_tmp = img_tmp(101:end-50,51:end-50);
        
        subplot2(2,3,1);
        imshow2(img_tmp,[0,10]);
        
        h = fspecial('average', [7 7]);
        img_tmp2 = filter2(h, img_tmp);
        
        subplot2(2,3,2);
        imshow2(img_tmp2,[0,10]); hold all;
        
        bwImg = imbinarize(img_tmp2,.9);
        skeletonImg = bwmorph(bwImg, 'skel', Inf);
        
        subplot2(2,3,3);
        imshow2(bwImg);        
        subplot2(2,3,4);
        imshow(skeletonImg);
        
        cleanedSkeleton = skeletonImg;
        % cleanedSkeleton = bwmorph(cleanedSkeleton, 'spur', 3);  % Adjust 5 for spurs of different lengths
        % branchPoints = bwmorph(cleanedSkeleton, 'branchpoints');
        % cleanedSkeleton = cleanedSkeleton & ~branchPoints;  % Remove branch points
        
        cc = bwconncomp(cleanedSkeleton);
        props = regionprops(cc, 'Area', 'ConvexHull', 'MajorAxisLength');
        rem = [props.Area] < 30;
        clear straightness;
        for i = 1:numel(props)
            skeletonLength = props(i).MajorAxisLength;            
            convexHull = props(i).ConvexHull;
            convexPerimeter = sum(sqrt(sum(diff([convexHull; convexHull(1,:)]).^2, 2)));
            straightness(i) = skeletonLength / convexPerimeter;
        end
        rem = find(rem | straightness<.5);
        for i = 1:length(rem)
            cleanedSkeleton(cc.PixelIdxList{rem(i)}) = 0;  % Set small components to 0
        end
        cc = bwconncomp(cleanedSkeleton);
        segmentedImg = labelmatrix(cc);
        subplot2(2,3,5);
        imshow(label2rgb(segmentedImg, 'jet', 'k', 'shuffle'));  hold on;
        title(num2str(cc.NumObjects));
                
        % merging
        props = regionprops(cc, 'Area', 'Orientation', 'PixelList', 'Centroid');
        proximityThreshold = 15;  % Maximum distance between components to be considered for merging
        orientationThreshold = 5;  % Maximum difference in orientation (degrees) for merging
        for i = 1:length(props)
            for j = i+1:length(props)
        
                pixelList1 = props(i).PixelList;
                pixelList2 = props(j).PixelList;
                minDistance = Inf;
                for p1 = 1:size(pixelList1, 1)
                    for p2 = 1:size(pixelList2, 1)
                        distance = sqrt((pixelList1(p1,1) - pixelList2(p2,1))^2 + ...
                                        (pixelList1(p1,2) - pixelList2(p2,2))^2);
                        if distance < minDistance
                            minDistance = distance;
                        end
                    end
                end
        
                centroid1 = props(i).Centroid;
                centroid2 = props(j).Centroid;               
                centroidOrientation = atan2d(-(centroid2(2) - centroid1(2)), centroid2(1) - centroid1(1));

                orientation1 = props(i).Orientation;
                orientation2 = props(j).Orientation;
                orientationDiff1 = abs(orientation1 - centroidOrientation);
                orientationDiff2 = abs(orientation2 - centroidOrientation);
                                     
                if minDistance < proximityThreshold && ...
                        orientationDiff1 < orientationThreshold && ...
                        orientationDiff2 < orientationThreshold
                    segmentedImg(segmentedImg == j) = i;
                end
            end
        end
        
        idx = setdiff(unique(segmentedImg),0);
        nline(exp,f) = numel(idx);
        
        % Display the final merged and labeled image
        subplot2(2,3,6);
        imshow(label2rgb(segmentedImg, 'jet', 'k', 'shuffle'));  hold on;
        title(['Merged Regions (', num2str(nline(exp,f)), ' objects)']);
        
        % Overlay indices on the merged regions
        for i = 1:nline(exp,f)
            [I,J] = find(segmentedImg==idx(i));
            text(mean(J), mean(I), num2str(idx(i)), 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
        end
        % pause;

        sfigure(100+exp);
        subplot2(3,3,f);
        img_tmp = repmat(1-rescale(img_tmp,'InputMax',10,'InputMin',-5),[1,1,3]) + double(label2rgb(segmentedImg, 'jet', 'k', 'shuffle'));
        imshow2(img_tmp);
    end
end

%% 
figure(2); clf;
ydat = mean(nline,2);
eydat = std(nline,[],2);
bar(ydat); hold all;
errorbar(ydat,eydat);

save('analysis2','img_proc','nline');




%%
pth = 'raw2_low mag\';
finfo = {dir([pth,'\*no tau_0*.tif']), dir([pth,'\*with tau_0*.tif'])};
finfo{1}(2:2:end) = []; % focus move images
npxr = 512; npxc = npxr;
img = cell(2,1); nfile = zeros(2,1);
for exp = 1:2
    nfile(exp) = numel(finfo{exp});
    img{exp} = zeros(npxr,npxc,nfile(exp),2);
    for f = 1:nfile(exp)
        disp([exp,f]);
        % img_tmp = imread([pth,finfo{exp}(f).name]);
        mov = loadTifStack16([pth,finfo{exp}(f).name]);
        img{exp}(:,:,f,1) = mean(mov(1:512,256+(1:512),:),3);
        img{exp}(:,:,f,2) = mean(mov(513:end,256+(1:512),:),3);
    end
end
clear mov;

%%
maxfig(11); clf;
Imax_list = [300,250];
% f_list = [5,5];
% rrange = 350+(-99:100); crange = 246+(-99:100);
rrange = 1:512; crange = 1:512;
for exp = 1:2
    for ch = 1:2
        for f = 1:4
            subplot(4,4,8*(exp-1)+4*(ch-1)+f);
            img_tmp = img{exp}(rrange,crange,f,ch);
            img_tmp2 = rescale(img_tmp,'InputMax',Imax_list(ch),'InputMin',150);
            imshow2(img_tmp2);
        end
    end
end

% save('analysis_raw2');