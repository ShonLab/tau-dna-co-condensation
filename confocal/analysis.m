%%
% rtpth = 'unsynced\1. Prometaphase\'; fname_base = 'sample';
% rtpth = 'synced\1. Double Thymidine block - HEK293\1. Prometaphase\'; fname_base = 'Sample ';
% rtpth = 'synced\2. Double Thymidine block - SH-SY5Y\1. Prometaphase\'; fname_base = 'Sample ';
% rtpth = 'synced\1. Double Thymidine block - HEK293\2. Prometa-Metaphase\'; fname_base = 'Sample ';
% rtpth = 'synced\2. Double Thymidine block - SH-SY5Y\2. Prometa-Metaphase\'; fname_base = 'Sample ';
rtpth = 'synced\1. Double Thymidine block - HEK293\4. Anaphase\'; fname_base = 'Sample ';
% rtpth = 'synced\2. Double Thymidine block - SH-SY5Y\2. Prometa-Metaphase\'; fname_base = 'Sample ';
% load([rtpth,'analysis']);
%%
read = 1;
folderinfo = dir([rtpth,fname_base,'*']); nexp = numel(folderinfo);

if read
    finfo = cell(nexp,1);
    img = cell(nexp,1); nfile = zeros(nexp,1);
    % cnt = [];
end
maxfig(nexp:-1:1);
% clear ax;
for i = 1:nexp
% for i = 15
    disp(i);

    if read
        pth = [rtpth,folderinfo(i).name,'\'];
        finfo{i} = dir([pth,'*.tif']); nfile(i) = numel(finfo{i});
        [~,idx] = natsort({finfo{i}.name}); finfo{i} = finfo{i}(idx);
    end
    
    sfigure(i);
    for j = 1:nfile(i)
        if j == 1
            clf;
        end

        if read
            img{i} = cat(4,img{i},loadTifStack16([pth,finfo{i}(j).name]));
        end
        img_tmp = [];
        for ch = 1:3
            img_tmp = cat(3, img_tmp, rescale(img{i}(:,:,ch,j),'InputMin',0));
        end
        ch = 4;
        img_tmp = img_tmp + repmat(rescale(img{i}(:,:,ch,j),'InputMin',0),[1,1,3]);

        ax(i,j) = subplot2(2,5,j);
        imshow2(img_tmp); title([num2str(i),', ',finfo{i}(j).name(end-5:end-4)]); hold all;
        plot([0,250],[10,10],'w','linew',2);
    end

    [~,j_max] = max(squeeze(msum(img{i}(:,:,4,:),[1,2])));
    tmp = centroid(img{i}(:,:,4,j_max));
    
    % axes(ax(i,j_max));
    % plot(tmp(2),tmp(1),'yo','markersize',10);
    % [tmpx, tmpy] = getpts;
    % cnt = [cnt; [i*ones(numel(tmpx),1),tmpx,tmpy]];
end
% save([rtpth,'analysis'],'img','cnt','finfo','nfile');

%%
% load([rtpth,'analysis']);
ncel = size(cnt,1);
[img_crop,xpos,ypos] = deal(cell(ncel,1));
[nmol,j_max] = deal(zeros(ncel,1));
span = 200; span2 = 20;
% Imax_list = [3,3,1,5]*1e3;

PCC = zeros(ncel,6);
img_selected = zeros(2*span+1,2*span+1,4,ncel,'uint16');
img_centromere = zeros(2*span2+1,2*span2+1,4,ncel,2);
I_centromere = zeros(2*span2+1,4,ncel,2);

maxfig(101:102); clf(101:102);
[x_list,y_list] = meshgrid(1:2*span+1);
for cel = 1:ncel
% for cel = 1:3
    disp(cel);
       
    i = cnt(cel,1);
    rrange = round(cnt(cel,3))+(-span:span);
    crange = round(cnt(cel,2))+(-span:span);
    
    sfigure(cel);
    for j = 1:nfile(i)
        if j == 1
            clf;
        end

        for ch = 1:4
            img_tmp = double(img{i}(rrange,crange,ch,j));
            img_tmp_sub = img_tmp(121:end-120,121:end-120);
            Imax = median(img_tmp_sub(:)) + 4*std(img_tmp_sub(:));
            img_crop{cel}(:,:,ch,j) = rescale(img_tmp,'InputMax',Imax);
        end
        
        ax(cel,j) = subplot2(2,5,j);
        imshow2(img_crop{cel}(:,:,[1,2,4],j)); title([num2str(cel),', ',finfo{i}(j).name(end-5:end-4)]); hold all;
        plot([0,250],[10,10],'w','linew',2);
    end
    [~,j_max(cel)] = max(squeeze(msum(img{i}(rrange(121:end-120),crange(121:end-120),4,:),[1,2]))); % maximum intensity in the centromere image
    j_max(19) = 3;
    img_selected(:,:,:,cel) = img{i}(rrange,crange,:,j_max(cel));
    
    sfigure(101);
    nrow = floor(sqrt(ncel/2)); ncol = ceil(ncel/nrow);
    subplot2(nrow,ncol,cel);
    imshow2(img_crop{cel}(:,:,[1,2,4],j_max(cel))); hold all;
    idx = 1;
    for ch1 = 1:3
        img1 = rescale(img_selected(:,:,ch1,cel));
        for ch2 = ch1+1:4
            img2 = rescale(img_selected(:,:,ch2,cel));
            numerator = sum((img1(:) - mean(img1)) .* (img2(:) - mean(img2(:))));
            denominator = sqrt(sum((img1(:) - mean(img1)).^2) .* sum((img2(:) - mean(img2)).^2));
            PCC(cel,idx) = numerator / denominator;
            idx = idx+1;
        end
    end    
    
    % locate centromeres
    img_tmp = img_selected(:,:,4,cel); threshold = .3*max(img_tmp(:));
    [nmol(cel),xpos{cel},ypos{cel}] = countSM(img_tmp, threshold,0);
    sel = hypot(xpos{cel}-(span+1),ypos{cel}-(span+1)) < 130;
    xpos{cel} = xpos{cel}(sel);
    ypos{cel} = ypos{cel}(sel);
    nmol(cel) = numel(xpos{cel});    
    plot(xpos{cel},ypos{cel},'m.','markersize',5);
    title([num2str(cel),', ',finfo{i}(j_max(cel)).name(end-5:end-4),', ',num2str(nmol(cel))]); hold all;

    % crop images around the centromeres
    for n = 1:nmol(cel)
        % dist = hypot(x_list(:)-xpos{cel}(n),y_list(:)-ypos{cel}(n));
        r0 = round(ypos{cel}(n)); c0 = round(xpos{cel}(n));
        rrange2 = r0 + (-span2:span2); crange2 = c0 + (-span2:span2);
        rrange3 = r0 + (-2*span2:2*span2); crange3 = c0 + (-2*span2:2*span2);
        for ch = 1:4
            img_tmp = rescale(img_selected(:,:,ch,cel));
            img_tmp_sub = img_tmp(rrange2,crange2);
            img_centromere(:,:,ch,cel,1) = img_centromere(:,:,ch,cel,1) + img_tmp_sub;
            I_centromere(:,ch,cel,1) = I_centromere(:,ch,cel,1) + mmean(img_tmp_sub(span2+1+(-1:1),:),1)';

            % orientation
            if ch == 1
                cnt_tmp = centroid(img_tmp(r0+(-5:5),c0+(-5:5)))-6;
                % figure(201); clf;
                % imshow2(img_tmp(rrange3,crange3)); hold on;
                % plot(2*span2+1,2*span2+1,'r.');
                % plot(2*span2+1+cnt_tmp(2),2*span2+1+cnt_tmp(1),'g.');
                % title(num2str([cel,n]));
                % pause;
            end
            img_tmp = imrotate(img_tmp(rrange3,crange3),-atan2d(-cnt_tmp(2),cnt_tmp(1)),'bicubic','crop');
            img_tmp_sub = img_tmp(span2+1:end-span2,span2+1:end-span2);
            img_centromere(:,:,ch,cel,2) = img_centromere(:,:,ch,cel,2) + img_tmp_sub;
            I_centromere(:,ch,cel,2) = I_centromere(:,ch,cel,2) + mmean(img_tmp_sub(span2+1+(-1:1),:),1)'; 
        end
    end
    for ch = 1:4
        img_centromere(:,:,ch,cel,1) = rescale(img_centromere(:,:,ch,cel,1));
        img_centromere(:,:,ch,cel,2) = rescale(img_centromere(:,:,ch,cel,2));
    end

    sfigure(102);
    nrow = floor(sqrt(ncel)); ncol = ceil(ncel*2/nrow);
    subplot2(nrow,ncol,cel);
    imshow2(img_centromere(:,:,4,cel),[]);
    subplot2(nrow,ncol,cel+ncel);
    imshow2(img_centromere(:,:,1:3,cel));
    title([num2str(cel),', ',finfo{i}(j_max(cel)).name(end-5:end-4),', ',num2str(nmol(cel))]); hold all;
end

%%
maxfig(103); clf;
for ch = 1:4
    subplot(3,4,ch);
    imshow2(mean(img_centromere(:,:,ch,:,1),4),[]); hold all;
    subplot(3,4,ch+4);
    imshow2(mean(img_centromere(:,:,ch,:,2),4),[]); hold all;
end

subplot(3,6,13);
plotSpread(PCC);

subplot(3,6,14);
histogram(nmol,0:5:60);

for i = 1:2
    subplot(3,6,14+i);
    imshow2(mean(img_centromere(:,:,[1,2,4],:,i),4));
    
    subplot(3,6,16+i);
    % plot(-span2:span2,mean(I_centromere(:,:,:,i),3));
    xdat = -span2:span2;
    Idat_all = squeeze(mean(img_centromere(span2+1+(-1:1),:,:,:,i),1));
    Idat_avg = mean(Idat_all,3);
    Idat_err = std(Idat_all,[],3);
    plot(xdat,Idat_avg(:,[1,2,4])); hold on;
    % errorbar(xdat,Idat_avg(:,[1,2,4]),Idat_err(:,[1,2,4]));
end

%%
saveas(101,[rtpth,'images.fig']);
saveas(102,[rtpth,'centromere.fig']);
saveas(103,[rtpth,'results.fig']);
save([rtpth,'analysis2'],'img_selected','img_centromere','I_centromere','nmol','PCC');