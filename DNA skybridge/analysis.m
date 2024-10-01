clear all; close all;
%% adjustable
% experimental condition
framerate=0.1;  %second
scale=108*10^(-9);  %meter

%parameters for colocalization
target_frame=1000;
threshold=[200, 200]; max_off=2; disp_range=[4000,1000]; 

%parameters for delete skybridge
barwidth=12; period=16;


% arrow_amplifier=30;
sz=2;

len=30; %minimum length that the particle reside in ROI
maxdisp=3; %the maximum distance particles can move after one frame


r_folder='result'; %folder to save results
mkdir(r_folder)
%%
plist=dir;
isSubfolder=[plist.isdir]; isSubfolder(1:2)=0;
path={plist(isSubfolder).name};
path(contains(path, 'result'))=[];

for p=1:size(path,2)
    datainfo=dir(path{p});
    fname{p}={datainfo.name};
    fname{p}(~contains(fname{p},'.tif'))=[];
    nfile{p}=numel(fname{p});


    for i=1:nfile{p}
        co_D_data{p}{i}={}; D_data{p}{i}={}; co_alpha_data{p}{i}={}; alpha_data{p}{i}={}; co_ccf{p}{i}={}; ccf{p}{i}={};

       img=strcat(path{p},'\',fname{p}{i});
       info=imfinfo(img);

       %image loading
       for imageNumber=1:size(info,1)
           img_raw(:,:,imageNumber)=imread(img,'index',imageNumber);
       end

       bkg=medfilt2(img_raw(:,:,end),round(size(img_raw(:,:,1))/10));

       %colocalization
        if ~isfile('tform.mat') && p==1 & i==1
        [tform,comp_img]=find_tform(img_raw(:,:,1)-bkg,disp_range);
        rname=strcat('tform result.tif'); delete(rname)
        writeTIFrgb(comp_img,rname);
        save('tform.mat', 'tform')
        else load('tform.mat') 
        end

        %vertical line deletion
        mask1=uint16(vertical_line_detection(img_raw(1+end/2:end,:,end)-bkg(1+end/2:end,:),barwidth,period)); 
        mask2=imwarp(mask1, invert(tform), 'OutputView', imref2d(size(mask1)));
        mask=[mask2;mask1];

        bkg3d=repmat(bkg,1,1,size(img_raw,3)); mask3d=repmat(mask,1,1,size(img_raw,3));
        img_mask=img_raw; img_mask(~mask3d)=bkg3d(~mask3d);
          
       [up,down]=dualviewer_merger(img_mask,tform,strcat(r_folder,'\composite_',fname{p}{i}));
       
       co_pos{p}{i}=colocal_1d(strcat(r_folder,'\colocal_',fname{p}{i},'.png'),up, down,threshold,max_off,disp_range);
       
       co_tr=[];
       for imageNumber=1:size(info,1)
           co_tr=[co_tr; [co_pos{p}{i}{imageNumber} repmat(imageNumber,size(co_pos{p}{i}{imageNumber},1),1)]];
       end

       %tracking from upper channel
       tr=particle_tracking(up, threshold(1), maxdisp, len);
       if isempty(tr) continue; end

       [~,it,ic]=intersect(tr(:,1:3),co_tr,'rows');
       co_tr(ic,4)=tr(it,4); co_tr(co_tr(:,4)==0,:)=[];

       if ~isempty(co_tr)
           c_ids=unique(co_tr(:,4));
           for co_id=1:numel(c_ids)
               co_tr=union(co_tr,tr(tr(:,4)==c_ids(co_id),:),'rows');
           end
           if ~isempty(co_tr)
               co_D{p}{i}=CalD_1D(1,co_tr, framerate, scale,strcat(r_folder,'\colocal msd ',fname{p}{i})); 
           end
       else
           co_D{p}{i}={};
       end

       tr=setdiff(tr,co_tr,'rows');
       if ~isempty(co_tr)
           nonco_D{p}{i}=CalD_1D(1,tr, framerate, scale,strcat(r_folder,'\noncolocal msd ',fname{p}{i}));
       else
           nonco_D{p}{i}={};
       end

  
       if ~isempty(co_D{p}{i})
          c_ids=unique(co_tr(:,4));
          for id=1:numel(c_ids)
              roi=sortrows(co_tr(find(co_tr(:,4)==c_ids(id)),:),3);
              y=round(roi(1,2));

              if nnz(~mask1(y,round(max(roi(1,1)-5,1):min(roi(1,1)+5,size(up,2))))) continue; end
              co_D_data{p}{i}{id,1}=co_D{p}{i}(find(co_D{p}{i}(:,3)==c_ids(id)),1); 
              co_alpha_data{p}{i}{id,1}=co_D{p}{i}(find(co_D{p}{i}(:,3)==c_ids(id)),2);

              x=round(min(roi(:,1))-2:max(roi(:,1))+2);
              
              t1=roi(1,3); t2=roi(end,3);

              k1=kymo(up(y,x,1:t2));
              k2=kymo(down(y,x,1:t2));

              ccf_temp=corrcoef(double(k1(t1:t2,:)),double(k2(t1:t2,:)));
              co_ccf{p}{i}{id,1}=ccf_temp(2,1);          

              fig = figure('Visible', 'off');
              ax1 = subplot(121);
              imagesc(k1);
              ax1.YAxis.Limits(1)=ax1.YAxis.Limits(1)+t1-1;
              title('tau')

              ax2 = subplot(122);
              imagesc(k2); 
              ax2.YAxis.Limits(1)=ax2.YAxis.Limits(1)+t1-1;
              hold on;
              title('DNA')
              
              linkaxes([ax1,ax2],'xy');
              sgtitle([fname{p}{i}, newline 'id=' num2str(c_ids(id)), ', ccf=' num2str(ccf_temp(2,1)), newline 'logD=', num2str(log10(co_D_data{p}{i}{id,1}*10^12)), 'um^2/s, alpha=', num2str(co_alpha_data{p}{i}{id,1})]);

              saveas(fig,strcat(r_folder,'\kymo colocalized_p=',num2str(p),'_i=',num2str(i),'_id=', num2str(c_ids(id)), '.png'))
          end
       else
           co_D_data{p}{i}=[]; 
           co_alpha_data{p}{i}=[];
           co_ccf{p}{i}=[];
       end



      if ~isempty(nonco_D{p}{i})
          ids=unique(tr(:,4));
          for id=1:numel(ids)
              roi=sortrows(tr(find(tr(:,4)==ids(id)),:),3);
              y=round(roi(1,2));

              if nnz(~mask1(y,round(max(roi(1,1)-5,1):min(roi(1,1)+5,size(up,2))))) continue; end
              D_data{p}{i}{id,1}=nonco_D{p}{i}(find(nonco_D{p}{i}(:,3)==ids(id)),1); 
              alpha_data{p}{i}{id,1}=nonco_D{p}{i}(find(nonco_D{p}{i}(:,3)==ids(id)),2);

              x=round(min(roi(:,1))-2:max(roi(:,1))+2);
              
              t1=roi(1,3); t2=roi(end,3);

              k1=kymo(up(y,x,1:t2));
              k2=kymo(down(y,x,1:t2));

              ccf_temp=corrcoef(double(k1(t1:t2,:)),double(k2(t1:t2,:)));
              ccf{p}{i}{id,1}=ccf_temp(2,1);          

              fig = figure('Visible', 'off');
              ax1 = subplot(121);
              imagesc(k1);
              ax1.YAxis.Limits(1)=ax1.YAxis.Limits(1)+t1-1;
              title('tau')

              ax2 = subplot(122);
              imagesc(k2); 
              ax2.YAxis.Limits(1)=ax2.YAxis.Limits(1)+t1-1;
              hold on;
              title('DNA')
              
              linkaxes([ax1,ax2],'xy');
              sgtitle([fname{p}{i}, newline 'id=' num2str(ids(id)), ', ccf=' num2str(ccf_temp(2,1)), newline 'logD=', num2str(log10(D_data{p}{i}{id,1}*10^12)), 'um^2/s, alpha=', num2str(alpha_data{p}{i}{id,1})]);

              saveas(fig,strcat(r_folder,'\kymo_p=',num2str(p),'_i=',num2str(i),'_id=', num2str(c_ids(id)), '.png'))
          end
       else
           D_data{p}{i}=[]; 
           alpha_data{p}{i}=[];
           ccf{p}{i}=[];
      end

       n1=numel(co_D_data{p}{i}); n2=numel(D_data{p}{i});

       if n1+n2==0; co_ratio{p}{i}=0;
       else co_ratio{p}{i}=n1/(n1+n2); end


        traj_nonco=trajectory(size(up),tr); %non-colocal
        traj_co=trajectory(size(up),co_tr);  %colocalized
        
        delete([r_folder '\traj_' fname{p}{i}]);
        for imageNumber=1:size(info,1)
            comp_img=cat(3,up(:,:,imageNumber),traj_nonco(:,:,imageNumber),traj_co(:,:,imageNumber));
            writeTIFrgb(comp_img,[r_folder '\traj_' fname{p}{i}]);
        end


        clear filtered_img img_raw tr particles roi bkg3d bkg up down img_mask mask3d mask1 mask2 mask
    end
end

save('data.mat')
%% plot-colocal ratio
clear;
load('result_best\data.mat');

labels={'colocal ratio'};
clear data
figure('Position',[700 300 300 300]);
data=[];
for p=1:size(co_ccf,2)
    data=[data; cell2mat(co_ratio{p}')];
end

bar(mean(data)); hold on;
errorbar(mean(data),std(data))
scatter(rand([numel(data),1])+0.5,data,12,'k','filled'); hold off;

title({'colocalized tau / whole tau'})

%% plot-total diffusion coeff
clear data

colocalized=[]; non_colocal=[];
for p=1:size(co_ccf,2)-1
    for i=1:nfile{p}
        colocalized=[colocalized;[log10(cell2mat(co_D_data{p}{i})*10^12) cell2mat(co_alpha_data{p}{i}) cell2mat(co_ccf{p}{i})]];
        non_colocal=[non_colocal;[log10(cell2mat(D_data{p}{i})*10^12) cell2mat(alpha_data{p}{i}) cell2mat(ccf{p}{i})]];
    end
end

colocalized(find(colocalized(:,3)<0.4),:)=[];


data=[colocalized(:,1:2)];

figure;
h=scatterhist(data(:,1),data(:,2), 'MarkerSize',1); hold on;  

xbin=0.1; ybin=0.05;
hx = h(2).Children;
for i = 1:length(hx)
    if isa(hx(i), 'matlab.graphics.chart.primitive.Histogram')
        hx(i).Normalization = 'probability';
        hx(i).BinWidth=xbin;
        hx(i).DisplayStyle='bar';
    end
end

hy = h(3).Children;
for i = 1:length(hy)
    if isa(hy(i), 'matlab.graphics.chart.primitive.Histogram')
        hy(i).Normalization = 'probability';
        hy(i).BinWidth=ybin;
        hy(i).DisplayStyle='bar';
    end
end
h(2).YLim = [0, 0.1];
h(3).YLim = [0, 0.5];
h(2).Visible = 'on';
h(3).Visible = 'on';
ylabel(h(2),'probability');
ylabel(h(3),'probability');
ylabel('alpha')
xlabel('log D (μm^2/s)')
title([num2str(size(colocalized,1)) ' colocalized spots, '])

figure;
bin=0.2;
histogram(colocalized(:,1),'Binwidth',bin,'Normalization','probability'); hold on
histogram(non_colocal(:,1),'Binwidth',bin,'Normalization','probability'); hold off
ylabel('probability')
xlabel('log D (um^2/s)');
title([num2str(size(colocalized,1)) ' colocalized spots, ' num2str(size(non_colocal,1)) ' non-colocalized spots'])
legend({'colocalized','noncolocalized'})

%%
clear data
labels={'colocalized','noncolocalized'};
xtick=1:2;

figure('Position',[700 300 200 300]);

data={colocalized(:,2), non_colocal(:,2)};

for p=1:2
    tmpdat=data{p}(:,1);
    bar(xtick(p),mean(tmpdat),0.5); hold on;
    errorbar(xtick(p),mean(tmpdat),std(tmpdat))
    scatter(ones(size(tmpdat))*xtick(p),tmpdat,12,'k','filled'); hold on;
end

xticks(xtick);
xticklabels(labels);
ylabel('alpha')

hold off

%%
clear data
labels={'colocalized','noncolocalized'};
xtick=1:2;

figure('Position',[700 300 200 300]);

data={colocalized(:,3), non_colocal(:,3)};

for p=1:2
    tmpdat=data{p}(:,1);
    bar(xtick(p),mean(tmpdat),0.3); hold on;
    errorbar(xtick(p),mean(tmpdat),std(tmpdat))
    scatter(ones(size(tmpdat))*xtick(p),tmpdat,12,'ko'); hold on;
end

hold off

xticks(xtick);
xticklabels(labels);
ylabel('correlation coefficients')







