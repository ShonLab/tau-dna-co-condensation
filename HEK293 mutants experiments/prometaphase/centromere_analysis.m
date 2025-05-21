%%
% rtpth = 'unsynced\1. Prometaphase\'; fname_base = 'sample';
% rtpth = 'synced\1. Double Thymidine block - HEK293\1. Prometaphase\'; fname_base = 'Sample ';
% rtpth = 'synced\2. Double Thymidine block - SH-SY5Y\1. Prometaphase\'; fname_base = 'Sample ';
% rtpth = 'synced\1. Double Thymidine block - HEK293\2. Prometa-Metaphase\'; fname_base = 'Sample ';
% rtpth = 'synced\2. Double Thymidine block - SH-SY5Y\2. Prometa-Metaphase\'; fname_base = 'Sample ';
% rtpth = 'S262D\'; fname_base = 'Sample';
rtpth = 'Double\'; fname_base = 'Sample';
% rtpth = 'synced\2. Double Thymidine block - SH-SY5Y\2. Prometa-Metaphase\'; fname_base = 'Sample ';
% load([rtpth,'analysis']);

%%
% read = 1;
% folderinfo = dir([rtpth,fname_base,'*']); nexp = numel(folderinfo);
% 
% if read
%     finfo = cell(nexp,1);
%     img = cell(nexp,1); nfile = zeros(nexp,1);
%     cnt = [];
% end
% 
% %
% 
% maxfig(nexp:-1:1);
% % clear ax;
% for i = 1:nexp
% % for i = 15
%     disp(i);
% 
%     if read
%         pth = [rtpth,folderinfo(i).name,'\'];
%         finfo{i} = dir([pth,'*.tif']); nfile(i) = numel(finfo{i});
%         [~,idx] = natsort({finfo{i}.name}); finfo{i} = finfo{i}(idx);
%     end
% 
%     sfigure(i);
%     for j = 1:nfile(i)
%         if j == 1
%             clf;
%         end
% 
%         if read
%             img{i} = cat(4,img{i},loadTifStack16([pth,finfo{i}(j).name]));
%         end
%         img_tmp = [];
%         for ch = [1,2,4]
%             img_tmp = cat(3, img_tmp, rescale(img{i}(:,:,ch,j),'InputMin',0));
%         end
%         ch = 3;
%         img_tmp = img_tmp + repmat(rescale(img{i}(:,:,ch,j),'InputMin',0),[1,1,3]);
% 
%         ax(i,j) = subplot2(1,1,j);
%         imshow2(img_tmp); title([num2str(i),', ',finfo{i}(j).name(end-5:end-4)]); hold all;
%         plot([0,250],[10,10],'w','linew',2);
%     end
% 
%     [~,j_max] = max(squeeze(msum(img{i}(:,:,4,:),[1,2])));
%     tmp = centroid(img{i}(:,:,4,j_max));
% %     
%     axes(ax(i,j_max));
%     plot(tmp(2),tmp(1),'yo','markersize',10);
%     [tmpx, tmpy] = getpts;
%     cnt = [cnt; [i*ones(numel(tmpx),1),tmpx,tmpy]];
% end

%%
% save([rtpth,'double_analysis'],'img','cnt','finfo','nfile');

%%
% load([rtpth,'s262d_analysis']);
load([rtpth,'double_analysis']);
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
% for cel = [1,2,3,4,5,8] % S262D
% for cel = [2,3,4,7,8] % S262D
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
        imshow2(img_crop{cel}(:,:,[1,2,3],j)); title([num2str(cel),', ',finfo{i}(j).name(end-5:end-4)]); hold all;
        plot([0,250],[10,10],'w','linew',2);
    end
    [~,j_max(cel)] = max(squeeze(msum(img{i}(rrange(121:end-120),crange(121:end-120),3,:),[1,2]))); % maximum intensity in the centromere image
    j_max(19) = 3;
    img_selected(:,:,:,cel) = img{i}(rrange,crange,:,j_max(cel));
    
    sfigure(101);
    nrow = floor(sqrt(ncel/2)); ncol = ceil(ncel/nrow);
    subplot2(nrow,ncol,cel);
    imshow2(img_crop{cel}(:,:,[1,2,3],j_max(cel))); hold all;
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
    img_tmp = img_selected(:,:,3,cel); threshold = .3*max(img_tmp(:));
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
    imshow2(img_centromere(:,:,3,cel),[]);
    subplot2(nrow,ncol,cel+ncel);
    imshow2(img_centromere(:,:,[1,2,4],cel));
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
    imshow2(mean(img_centromere(:,:,[1,2,3],:,i),4));
    
    subplot(3,6,16+i);
    % plot(-span2:span2,mean(I_centromere(:,:,:,i),3));
    xdat = -span2:span2;
    Idat_all = squeeze(mean(img_centromere(span2+1+(-1:1),:,:,:,i),1));
    Idat_avg = mean(Idat_all,3);
    Idat_err = std(Idat_all,[],3);
    plot(xdat,Idat_avg(:,[1,2,3])); hold on;
    % errorbar(xdat,Idat_avg(:,[1,2,4]),Idat_err(:,[1,2,4]));
end

%%
% saveas(101,[rtpth,'images.fig']);
% saveas(102,[rtpth,'centromere.fig']);
% saveas(103,[rtpth,'results.fig']);
% save([rtpth,'analysis2'],'img_selected','img_centromere','I_centromere','nmol','PCC');

%%
figure(104)
% Pearson correlation calculation for specific channel combinations in I_centromere
PCC_crop = zeros(ncel, 2); % Initialize Pearson Correlation Coefficients for two channel combinations (ch1,2 and ch2,4)

% Loop over each cell
for cel = 1:ncel
    % Calculate Pearson correlation between channels 1 and 2
    img1 = I_centromere(:,1,cel,2); % Channel 1 data
    img2 = I_centromere(:,2,cel,2); % Channel 2 data
    numerator = sum((img1 - mean(img1)) .* (img2 - mean(img2)));
    denominator = sqrt(sum((img1 - mean(img1)).^2) * sum((img2 - mean(img2)).^2));
    PCC_crop(cel,1) = numerator / denominator; % Store the result for channels 1 and 2

    % Calculate Pearson correlation between channels 2 and 4
    img1 = I_centromere(:,2,cel,2); % Reuse channel 2 data
    img2 = I_centromere(:,3,cel,2); % Channel 4 data
    numerator = sum((img1 - mean(img1)) .* (img2 - mean(img2)));
    denominator = sqrt(sum((img1 - mean(img1)).^2) * sum((img2 - mean(img2)).^2));
    PCC_crop(cel,2) = numerator / denominator; % Store the result for channels 2 and 4
end

%
sampleNames={'Tubulin-Tau','Tau-Centromere','Cropped Tubulin-Tau','Cropped Tau-Centromere'};
% Display the Pearson Correlation Coefficients

plotSpread([PCC(:,[1,5]),PCC_crop(:,:)]);
title('PCC',' ')
% x축 레이블 설정
set(gca, 'XTickLabel', sampleNames);




%% Calculate PCC from centromere region
% 데이터 로드
figure(105)
for cel = 1:ncel
    img1 = img_centromere(:,:,1,cel,2); % Channel 1 data
    img2 = img_centromere(:,:,2,cel,2); % Channel 2 data
    
    % Pearson correlation 계산
    R = corrcoef(img1(:), img2(:)); % 행렬을 벡터로 변환하여 계산
    PCC_value = R(1,2); % R은 2x2 행렬이며, 대각선 이외의 값이 correlation coefficient
    
    % 결과 출력
%     fprintf('Pearson Correlation Coefficient: %f\n', PCC_value);
    
    % 산점도 그리기
    subplot(3,6,cel)
    scatter(img1(:), img2(:), 1, 'filled'); % '10'은 마커의 크기를 나타냄
    title(sprintf('PCC of Tubulin-Tau (PCC = %.2f)', PCC_value));
    xlabel('Channel 1 Pixel Values');
    ylabel('Channel 2 Pixel Values');
    grid on;
    
    % 추가적인 선형 피팅 라인 그리기
    hold on;
    fitLine = polyfit(img1(:), img2(:), 1); % 1차 다항식으로 데이터에 적합
    refLine = polyval(fitLine, [min(img1(:)) max(img1(:))]);
    plot([min(img1(:)) max(img1(:))], refLine, 'r--', 'LineWidth', 2);
    hold off;
    PCC_pixel_crop(cel,1)=PCC_value;
end

figure(106)
for cel = 1:ncel
    img1 = img_centromere(:,:,2,cel,2); % Channel 1 data
    img2 = img_centromere(:,:,3,cel,2); % Channel 2 data
    
    % Pearson correlation 계산
    R = corrcoef(img1(:), img2(:)); % 행렬을 벡터로 변환하여 계산
    PCC_value = R(1,2); % R은 2x2 행렬이며, 대각선 이외의 값이 correlation coefficient
    
    % 결과 출력
%     fprintf('Pearson Correlation Coefficient: %f\n', PCC_value);
    
    % 산점도 그리기
    subplot(3,6,cel)
    scatter(img1(:), img2(:), 1, 'filled'); % '10'은 마커의 크기를 나타냄
    title(sprintf('PCC of Tau-Cent. (PCC = %.2f)', PCC_value));
    xlabel('Channel 1 Pixel Values');
    ylabel('Channel 2 Pixel Values');
    grid on;
    
    % 추가적인 선형 피팅 라인 그리기
    hold on;
    fitLine = polyfit(img1(:), img2(:), 1); % 1차 다항식으로 데이터에 적합
    refLine = polyval(fitLine, [min(img1(:)) max(img1(:))]);
    plot([min(img1(:)) max(img1(:))], refLine, 'r--', 'LineWidth', 2);
    hold off;
    PCC_pixel_crop(cel,2)=PCC_value;
end


%% Calculate PCC from selected region
% 데이터 로드
figure(107)
for cel = 1:ncel
    img1 = img_selected(:,:,1,cel); % Channel 1 data
    img2 = img_selected(:,:,2,cel); % Channel 2 data
    img1 = rescale(double(img1));
    img2 = rescale(double(img2));
    
    % Pearson correlation 계산
    R = corrcoef(img1(:), img2(:)); % 행렬을 벡터로 변환하여 계산
    PCC_value = R(1,2); % R은 2x2 행렬이며, 대각선 이외의 값이 correlation coefficient
    
    % 결과 출력
%     fprintf('Pearson Correlation Coefficient: %f\n', PCC_value);
    
    % 산점도 그리기
    subplot(3,6,cel)
    scatter(img1(:), img2(:), 1, 'filled'); % '10'은 마커의 크기를 나타냄
    title(sprintf('PCC of Tubulin-Tau (PCC = %.2f)', PCC_value));
    xlabel('Channel 1 Pixel Values');
    ylabel('Channel 2 Pixel Values');
    grid on;
    
    % 추가적인 선형 피팅 라인 그리기
    hold on;
    fitLine = polyfit(img1(:), img2(:), 1); % 1차 다항식으로 데이터에 적합
    refLine = polyval(fitLine, [min(img1(:)) max(img1(:))]);
    plot([min(img1(:)) max(img1(:))], refLine, 'r--', 'LineWidth', 2);
    hold off;
    PCC_pixel_whole(cel,1)=PCC_value;
end

figure(108)
for cel = 1:ncel
    img1 = img_selected(:,:,2,cel); % Channel 1 data
    img2 = img_selected(:,:,3,cel); % Channel 2 data
    img1 = rescale(double(img1));
    img2 = rescale(double(img2));
    
    % Pearson correlation 계산
    R = corrcoef(img1(:), img2(:)); % 행렬을 벡터로 변환하여 계산
    PCC_value = R(1,2); % R은 2x2 행렬이며, 대각선 이외의 값이 correlation coefficient
    
    % 결과 출력
%     fprintf('Pearson Correlation Coefficient: %f\n', PCC_value);
    
    % 산점도 그리기
    subplot(3,6,cel)
    scatter(img1(:), img2(:), 1, 'filled'); % '10'은 마커의 크기를 나타냄
    title(sprintf('PCC of Tau-Cent. (PCC = %.2f)', PCC_value));
    xlabel('Channel 1 Pixel Values');
    ylabel('Channel 2 Pixel Values');
    grid on;
    
    % 추가적인 선형 피팅 라인 그리기
    hold on;
    fitLine = polyfit(img1(:), img2(:), 1); % 1차 다항식으로 데이터에 적합
    refLine = polyval(fitLine, [min(img1(:)) max(img1(:))]);
    plot([min(img1(:)) max(img1(:))], refLine, 'r--', 'LineWidth', 2);
    hold off;
    PCC_pixel_whole(cel,2)=PCC_value;
end


%%
ff = figure(109);
ff.Position = [100 550 800 400];
% Pearson correlation calculation for specific channel combinations in I_centromere
% PCC_crop = zeros(ncel, 2); % Initialize Pearson Correlation Coefficients for two channel combinations (ch1,2 and ch2,4)


sampleNames={'Selected row T-T','Selected row T-C','Selected all T-T','Selected all T-C','Crop row T-T','Crop row T-C','Crop all T-T','Crop all T-C'};
% Display the Pearson Correlation Coefficients
subplot(1,2,1)
plotSpread([PCC(:,[1,5]),PCC_pixel_whole,PCC_crop,PCC_pixel_crop]);
title('PCC',' ')
% x축 레이블 설정
set(gca, 'XTickLabel', sampleNames);

sampleNames={'Crop row T-T','Crop row T-C','Crop all T-T','Crop all T-C'};
subplot(1,2,2)
plotSpread([PCC_crop,PCC_pixel_crop]);
title('PCC',' ')
ylim([-0.2 1])
% x축 레이블 설정
set(gca, 'XTickLabel', sampleNames);


%%

% 데이터 변수 준비
PCC_matrix = [PCC_crop,PCC_pixel_crop];

% 그룹 이름 설정
groupNames = {'row-avg Tau-Tubulin','row-avg Tau-Centromere','raw Tau-Tubulin','raw Tau-Centromere'};

% 그래프 설정
ff=figure(110);
ff.Position = [100 350 400 500];

% plotSpread 사용
plotSpread(PCC_matrix,'distributionColors', 'k', 'markerSize', 7);

% 데이터 행렬과 그룹 벡터 준비 (boxplot 용)
groupVector = repmat(1:size(PCC_matrix, 2), size(PCC_matrix, 1), 1);
groupVector = groupVector(:);
dataMatrix = PCC_matrix(:);
ylim([-1 1])
hold on; % 현재 축을 유지하고 추가 그래프 그리기
% boxplot 추가
h = boxplot(dataMatrix, groupVector, 'positions', 1:size(PCC_matrix, 2), 'Widths', 0.5, 'Colors', 'k', 'Symbol', 'k+');
lines = findobj(h, 'type', 'line'); % 박스플롯의 모든 선분(Line) 객체 찾기
set(lines, 'LineWidth', 1.5); % 모든 선분의 두께를 2로 설정
set(findobj(h, 'tag', 'Median'), 'Color', 'red', 'LineWidth', 2);

% 축 설정과 레이블
set(gca, 'XTickLabel', groupNames);
ylabel('Pearson Correlation Coefficient (PCC)');
title('PCC at centromere region');
hold off; % 추가 그리기 종료


%%
% 데이터 변수 준비
PCC_matrix = [PCC(:,[1,5]),PCC_pixel_whole];

% 그룹 이름 설정
groupNames = {'row-avg Tau-Tubulin','row-avg Tau-Centromere','raw Tau-Tubulin','raw Tau-Centromere'};

% 그래프 설정
ff=figure(111);
ff.Position = [500 350 400 500];
ax1 = subplot(1,1,1);

% plotSpread 사용
plotSpread(PCC_matrix,'distributionColors', 'k', 'markerSize', 7);
%
% 데이터 행렬과 그룹 벡터 준비 (boxplot 용)
groupVector = repmat(1:size(PCC_matrix, 2), size(PCC_matrix, 1), 1);
groupVector = groupVector(:);
dataMatrix = PCC_matrix(:);

hold on; % 현재 축을 유지하고 추가 그래프 그리기
% boxplot 추가
h = boxplot(dataMatrix, groupVector, 'positions', 1:size(PCC_matrix, 2), 'Widths', 0.5, 'Colors', 'k', 'Symbol', 'k+');
lines = findobj(h, 'type', 'line'); % 박스플롯의 모든 선분(Line) 객체 찾기
set(lines, 'LineWidth', 1.5); % 모든 선분의 두께를 2로 설정
set(findobj(h, 'tag', 'Median'), 'Color', 'red', 'LineWidth', 2);

% 축 설정과 레이블
set(gca, 'XTickLabel', groupNames);
ylabel('Pearson Correlation Coefficient (PCC)');
title('PCC at selected region');
hold off; % 추가 그리기 종료

sum(nmol)
mean(PCC_pixel_crop)
median(PCC_pixel_crop)

%%
ff=figure(115);
ff.Position = [500 350 350 500];
ch_list = [3,2,1];
color_list = {'b','g','r'};
ch_name = {'Centromere','Tau','Tubulin'};
for i = 1:3
    subplot(3,2,2*i-1);
    img = mean(img_centromere(:,:,ch_list(i),:,2), 4);
    
    % RGB 이미지 생성
    RGB = zeros(size(img, 1), size(img, 2), 3); % RGB 이미지 틀 생성
    switch i
        case 1
            RGB(:,:,3) = img; % Blue 채널에 할당
        case 2
            RGB(:,:,2) = img; % Green 채널에 할당
        case 3
            RGB(:,:,1) = img; % Red 채널에 할당
    end
    
    % 색상 이미지 표시
    imshow(RGB, []);

    % 선 그리기
    line([1, size(img, 2)], [size(img, 1)/2, size(img, 1)/2], 'Color', 0.9*[1 1 1], 'LineStyle', '--', 'LineWidth', 0.5); % 더 얇은 선
    
    % 이미지 크기에 따라 스케일 바 위치와 크기 조절
    imgSize = size(img);
    scaleBarLength = 11.628; % 이미지 너비의 20%
    scaleBarHeight = imgSize(1) * 0.05; % 이미지 높이의 5%
    scaleBarX = imgSize(2) * 0.1; % 왼쪽으로부터 10% 위치
    scaleBarY = imgSize(1) * 0.9; % 아래로부터 5% 위치
    
    % 스케일 바 그리기
    rectangle('Position', [scaleBarX, scaleBarY, scaleBarLength, scaleBarHeight], ...
              'FaceColor', 'w', 'EdgeColor', 'none');
    
    % 스케일 바 레이블
    text(scaleBarX + scaleBarLength / 1.5, scaleBarY - scaleBarHeight*2, ...
         '500 nm', 'Color', 'w', 'FontSize', 9, 'HorizontalAlignment', 'center');
    hold on;
    text(scaleBarX/1.5, scaleBarHeight*2.5, ...
         ch_name{i}, 'Color', 'w', 'FontSize', 10, 'HorizontalAlignment', 'left');
    hold off;
    
    subplot(3,2,2*i);
    % plot(-span2:span2,mean(I_centromere(:,:,:,i),3));
    xdat = -span2:span2;
    Idat_all = squeeze(mean(img_centromere(span2+1+(-1:1),:,:,:,2),1));
    Idat_avg = mean(Idat_all,3);
    Idat_err = std(Idat_all,[],3);
    plot(xdat*0.043,Idat_avg(:,ch_list(i)),color_list{i}); hold on;
    hold off;
    xlim([-0.75 0.75])
    ylim([0 1])
    ylabel('Intensity (a.u.)')
    
    if i==3
        xlabel('x position (um)')
    end
    % errorbar(xdat,Idat_avg(:,[1,2,4]),Idat_err(:,[1,2,4]));
end

%%
PCC_matrix = [PCC_pixel_whole(:,1),PCC_crop(:,2)];
groupVector = repmat(1:size(PCC_matrix, 2), size(PCC_matrix, 1), 1);
groupVector = groupVector(:);
dataMatrix = PCC_matrix(:);
groupNames = {'Tau-Tubulin','Tau-Centromere'};



ff = figure(121);
ff.Position = [100, 100, 700, 875];  % 전체 그림 크기 조정

% plotSpread와 Boxplot 결합
subplot(4, 3, [1,4]);
plotSpread(PCC_matrix,'distributionColors', 'k', 'markerSize', 7);
hold on;
boxplot(dataMatrix, groupVector, 'positions', 1:size(PCC_matrix, 2), 'Widths', 0.5, 'Colors', 'k', 'Symbol', 'k+');
lines = findobj(gca, 'type', 'line');
set(lines, 'LineWidth', 1.5);
set(findobj(gca, 'tag', 'Median'), 'Color', 'red', 'LineWidth', 2);
ylabel('Pearson Correlation Coefficient (PCC)');
set(gca, 'XTickLabel', groupNames);
% title('PCC at centromere region');
ylim([-0.2 1])
hold off;

% 각 채널별 이미지 및 인텐시티 플롯
for i = 1:3
    subplot(4,3,3*i);
    img = mean(img_centromere(:,:,ch_list(i),:,2), 4);
    
    % RGB 이미지 생성
    RGB = zeros(size(img, 1), size(img, 2), 3); % RGB 이미지 틀 생성
    switch i
        case 1
            RGB(:,:,3) = img; % Blue 채널에 할당
        case 2
            RGB(:,:,2) = img; % Green 채널에 할당
        case 3
            RGB(:,:,1) = img; % Red 채널에 할당
    end
    
    % 색상 이미지 표시
    imshow(RGB, []);

    % 선 그리기
    line([1, size(img, 2)], [size(img, 1)/2, size(img, 1)/2], 'Color', 0.9*[1 1 1], 'LineStyle', '--', 'LineWidth', 0.5); % 더 얇은 선
    
    % 이미지 크기에 따라 스케일 바 위치와 크기 조절
    imgSize = size(img);
    scaleBarLength = 11.628; % 이미지 너비의 20%
    scaleBarHeight = imgSize(1) * 0.05; % 이미지 높이의 5%
    scaleBarX = imgSize(2) * 0.1; % 왼쪽으로부터 10% 위치
    scaleBarY = imgSize(1) * 0.9; % 아래로부터 5% 위치
    
    % 스케일 바 그리기
    rectangle('Position', [scaleBarX, scaleBarY, scaleBarLength, scaleBarHeight], ...
              'FaceColor', 'w', 'EdgeColor', 'none');
    
    % 스케일 바 레이블
    text(scaleBarX + scaleBarLength /2, scaleBarY - scaleBarHeight*1.5, ...
         '500 nm', 'Color', 'w', 'FontSize', 11, 'HorizontalAlignment', 'center');
    hold on;
    text(scaleBarX/1.5, scaleBarHeight*2, ...
         ch_name{i}, 'Color', 'w', 'FontSize', 12, 'HorizontalAlignment', 'left');
    hold off;
    
    subplot(4,3,3*i-1);
    % plot(-span2:span2,mean(I_centromere(:,:,:,i),3));
    xdat = -span2:span2;
    Idat_all = squeeze(mean(img_centromere(span2+1+(-1:1),:,:,:,2),1));
    Idat_avg = mean(Idat_all,3);
    Idat_err = std(Idat_all,[],3);
    plot(xdat*0.043,Idat_avg(:,ch_list(i)),color_list{i}); hold on;
    hold off;
    xlim([-0.75 0.75])
    ylim([0 1])
    ylabel('Intensity (a.u.)')
    xlabel('x position (µm)')
    
end

% % plotSpread와 Boxplot 결합
% subplot(3, 3, [1 4 7]);
% plotSpread(PCC_matrix,'distributionColors', 'k', 'markerSize', 7);
% hold on;
% boxplot(dataMatrix, groupVector, 'positions', 1:size(PCC_matrix, 2), 'Widths', 0.5, 'Colors', 'k', 'Symbol', 'k+');
% lines = findobj(gca, 'type', 'line');
% set(lines, 'LineWidth', 1.5);
% set(findobj(gca, 'tag', 'Median'), 'Color', 'red', 'LineWidth', 2);
% ylabel('Pearson Correlation Coefficient (PCC)');
% set(gca, 'XTickLabel', groupNames);
% title('PCC at centromere region');
% ylim([-0.2 1])
% hold off;

% for i = 1:3
%     subplot(3,3,3*i-1);
%     img = mean(img_centromere(:,:,ch_list(i),:,2), 4);
%     
%     % RGB 이미지 생성
%     RGB = zeros(size(img, 1), size(img, 2), 3); % RGB 이미지 틀 생성
%     switch i
%         case 1
%             RGB(:,:,3) = img; % Blue 채널에 할당
%         case 2
%             RGB(:,:,2) = img; % Green 채널에 할당
%         case 3
%             RGB(:,:,1) = img; % Red 채널에 할당
%     end
%     
%     % 색상 이미지 표시
%     imshow(RGB, []);
% 
%     % 선 그리기
%     line([1, size(img, 2)], [size(img, 1)/2, size(img, 1)/2], 'Color', 0.9*[1 1 1], 'LineStyle', '--', 'LineWidth', 0.5); % 더 얇은 선
%     
%     % 이미지 크기에 따라 스케일 바 위치와 크기 조절
%     imgSize = size(img);
%     scaleBarLength = 11.628; % 이미지 너비의 20%
%     scaleBarHeight = imgSize(1) * 0.05; % 이미지 높이의 5%
%     scaleBarX = imgSize(2) * 0.1; % 왼쪽으로부터 10% 위치
%     scaleBarY = imgSize(1) * 0.9; % 아래로부터 5% 위치
%     
%     % 스케일 바 그리기
%     rectangle('Position', [scaleBarX, scaleBarY, scaleBarLength, scaleBarHeight], ...
%               'FaceColor', 'w', 'EdgeColor', 'none');
%     
%     % 스케일 바 레이블
%     text(scaleBarX + scaleBarLength /2, scaleBarY - scaleBarHeight*1.5, ...
%          '500 nm', 'Color', 'w', 'FontSize', 11, 'HorizontalAlignment', 'center');
%     hold on;
%     text(scaleBarX/1.5, scaleBarHeight*2, ...
%          ch_name{i}, 'Color', 'w', 'FontSize', 12, 'HorizontalAlignment', 'left');
%     hold off;
%     
%     subplot(3,3,3*i);
%     % plot(-span2:span2,mean(I_centromere(:,:,:,i),3));
%     xdat = -span2:span2;
%     Idat_all = squeeze(mean(img_centromere(span2+1+(-1:1),:,:,:,2),1));
%     Idat_avg = mean(Idat_all,3);
%     Idat_err = std(Idat_all,[],3);
%     plot(xdat*0.043,Idat_avg(:,ch_list(i)),color_list{i}); hold on;
%     hold off;
%     xlim([-0.75 0.75])
%     ylim([0 1])
%     ylabel('Intensity (a.u.)')
%     xlabel('x position (µm)')
%     
% end

%%


ff = figure(122);
ff.Position = [100, 100, 500, 500];  % 전체 그림 크기 조정

% 각 채널별 이미지 및 인텐시티 플롯
for i = 1:2
    subplot(2,2,2*i);
    img = mean(img_centromere(:,:,ch_list(i),:,2), 4);
    
    % RGB 이미지 생성
    RGB = zeros(size(img, 1), size(img, 2), 3); % RGB 이미지 틀 생성
    switch i
        case 1
            RGB(:,:,3) = img; % Blue 채널에 할당
        case 2
            RGB(:,:,2) = img; % Green 채널에 할당
        case 3
            RGB(:,:,1) = img; % Red 채널에 할당
    end
    
    % 색상 이미지 표시
    imshow(RGB, []);

    % 선 그리기
    line([1, size(img, 2)], [size(img, 1)/2, size(img, 1)/2], 'Color', 0.9*[1 1 1], 'LineStyle', '--', 'LineWidth', 0.5); % 더 얇은 선
    
    % 이미지 크기에 따라 스케일 바 위치와 크기 조절
    imgSize = size(img);
    scaleBarLength = 11.628; % 이미지 너비의 20%
    scaleBarHeight = imgSize(1) * 0.05; % 이미지 높이의 5%
    scaleBarX = imgSize(2) * 0.1; % 왼쪽으로부터 10% 위치
    scaleBarY = imgSize(1) * 0.9; % 아래로부터 5% 위치
    
    % 스케일 바 그리기
    rectangle('Position', [scaleBarX, scaleBarY, scaleBarLength, scaleBarHeight], ...
              'FaceColor', 'w', 'EdgeColor', 'none');
    
    % 스케일 바 레이블
    text(scaleBarX + scaleBarLength /2, scaleBarY - scaleBarHeight*1.5, ...
         '500 nm', 'Color', 'w', 'FontSize', 11, 'HorizontalAlignment', 'center');
    hold on;
    text(scaleBarX/1.5, scaleBarHeight*2, ...
         ch_name{i}, 'Color', 'w', 'FontSize', 12, 'HorizontalAlignment', 'left');
    hold off;
    
    subplot(2,2,2*i-1);
    % plot(-span2:span2,mean(I_centromere(:,:,:,i),3));
    xdat = -span2:span2;
    Idat_all = squeeze(mean(img_centromere(span2+1+(-1:1),:,:,:,2),1));
    Idat_avg = mean(Idat_all,3);
    Idat_err = std(Idat_all,[],3);
    plot(xdat*0.043,Idat_avg(:,ch_list(i)),color_list{i}); hold on;
    hold off;
    xlim([-0.75 0.75])
    ylim([0 1])
    ylabel('Intensity (a.u.)')
    xlabel('x position (µm)')
    
end


