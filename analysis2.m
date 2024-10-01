%% DNA envelope
% load('analysis');
p_list = {1:4, 5:6};
title_list = {'100 mM NaCl','25 mM NaCl'};
legend_list = {{'0 nM','50 nM','500 nM','5000 nM'}, {'0 nM','50 nM'}};

maxfig(1); clf; clear colororder;
colororder{1} = zeros(4,3); colororder{1}(:,1) = linspace(0,1,4);
colororder{2} = colororder{1}(1:2,:);
[sel1,sel2] = deal(cell(2,1));
rem = {{[54], [], [2,4], [4]}, {[], [8]}};
% 54: something on the att point
% 2,4: duplicates
% 4: duplicates
% 8: tau img clipped

for d = 1:numel(p_list)
    set(gcf, 'defaultaxesColorOrder', colororder{d});
    [xdat,ydat,edat,sel1{d},sel2{d}] = deal(cell(numel(p_list{d}),1));
    for pi = 1:numel(p_list{d})
        p = p_list{d}(pi);
        xdat{pi} = cat(1,att_dist{p}{:});
        ydat{pi} = cellfun(@(x) x.A*2*pxsz, cat(1,env_fitres{p}{:}));
        edat{pi} = cellfun(@(x) x.rmse, cat(1,env_gof{p}{:}));
        sel1{d}{pi} = edat{pi}<.9; % remove data that didn't fit the sine curve well
        % further discard bad ones
        if ~isempty(rem{d}{pi})
            sel1{d}{pi}(rem{d}{pi}) = false;
        end
        sel2{d}{pi} = xdat{pi}>8 & xdat{pi}<11; % limit DNA length
    end

    subplot(2,3,3*(d-1)+1);
    for pi = 1:numel(p_list{d})
        plot(xdat{pi}(sel1{d}{pi}),ydat{pi}(sel1{d}{pi}),'.'); hold all;
    end
    xlim([2,14]); ylim([0,2]);
    yxlabel('Envelope width (μm)','End-to-end distance (μm)');
    title(title_list{d});
    legend(legend_list{d});        

    for pid = 1:2
        subplot(2,3,3*(d-1)+1+pid);
        if pid == 1
            ydat_all = cellfun(@(y,s) y(s), ydat,sel1{d}, 'unif',0);
        else
            ydat_all = cellfun(@(y,s1,s2) y(s1&s2), ydat,sel1{d},sel2{d}, 'unif',0);
        end
        ydat_avg = cellfun(@(y) mean(y), ydat_all);
        ydat_std = cellfun(@(y) std(y), ydat_all);

        for pi = 1:numel(p_list{d})
            bar(pi,ydat_avg(pi),'facea',.6,'linestyle','none'); hold all;
        end
        errorbar(ydat_avg,ydat_std,'k','linestyle','none','CapSize',10);
        h = plotSpread(ydat_all,'xValues',1:numel(p_list{d}),'spreadWidth',.5, 'distributionColors',colororder{d});
        ylim([0,2]);
        set(gca,'xticklabel',legend_list{d});
        yxlabel('Envelope width (μm)','Tau concentration (nM)');
    end
end
saveas(gcf,'DNA envelope.jpg');

%% tau vs DNA
maxfig(2); clf;
for d = 1:numel(p_list)
    set(gcf, 'defaultaxesColorOrder', colororder{d});

    for pid = 1:6
        subplot(2,6,6*(d-1)+pid);
        ydat = cell(numel(p_list{d}),1);
        for pi = 1:numel(p_list{d})
            p = p_list{d}(pi);
            if pid <= 2
                for q = 1:nsubpth(p)
                    ydat{pi} = [ydat{pi}; cellfun(@(I1,I2) [mean(I1),mean(I2)], I_DNA{p}{q},I_tau{p}{q},'unif',0)];
                end
                ydat{pi} = cat(1,ydat{pi}{:});
            elseif pid <= 4
                ydat{pi} = cat(1,RMSd{p}{:});
            else
                ydat{pi} = cat(1,pearson{p}{:});
            end
        end
        ydat_all = cellfun(@(y,s1,s2) y(s1&s2,mod(pid-1,2)+1), ydat,sel1{d},sel2{d}, 'unif',0);
        ydat_avg = cellfun(@(y) nanmean(y), ydat_all, 'unif',0); ydat_avg = cat(1,ydat_avg{:});
        ydat_std = cellfun(@(y) nanstd(y), ydat_all, 'unif',0); ydat_std = cat(1,ydat_std{:});
    
        for pi = 1:numel(p_list{d})
            bar(pi,ydat_avg(pi),'facea',.6,'linestyle','none'); hold all;
        end
        errorbar(ydat_avg,ydat_std,'k','linestyle','none','CapSize',10);
        h = plotSpread(ydat_all,'xValues',1:numel(p_list{d}),'spreadWidth',.5, 'distributionColors',colororder{d});
        set(gca,'xticklabel',legend_list{d});
        if pid <= 2
            yxlabel('Mean intensity','Tau concentration (nM)');
        elseif pid <= 4
            ylim([0,2]);
            % yxlabel('Coefficient of variation','Tau concentration (nM)');
            yxlabel('RMS roughness','Tau concentration (nM)');
        else
            ylim([-.7,1]);
            yxlabel('Pearson correlation','Tau concentration (nM)');
        end
        
        if pid == 1 || pid == 3
            title('DNA fluorescence');
        elseif pid == 2 || pid == 4
            title('tau fluorescence');
        elseif pid == 5
            title('Forward correlation');
        else
            title('Reverse correlation');
        end
    end
end
saveas(gcf,'tau vs DNA.jpg');

%% sample images
% % find representative ones
% p=2; tmp = [cat(1,att_dist{p}{:}), cellfun(@(x) x.A*2*pxsz, cat(1,env_fitres{p}{:})), cat(1,pearson{p}{:})];
% d=1; pi=2; [find(sel1{d}{pi}&sel2{d}{pi}),tmp(sel1{d}{pi}&sel2{d}{pi},:)]

hex_DNA = '#7FFF00';
hex_tau = '#CC79A7';
cmap_DNA = hex2map(hex_DNA);
cmap_tau = hex2map(hex_tau);

q_list = {[1,1,1,1],[1,3]};
m_list = {[33,10,19,3],[49,9]};

outputVideo = VideoWriter('double-tether condensation.mp4', 'MPEG-4');
outputVideo.FrameRate = 10;
open(outputVideo);

figure(3); clf; set(gcf,'color',.1*ones(1,3));
for f = 1:50
    for d = 1:numel(p_list)
        for pi = numel(p_list{d}):-1:1
            p = p_list{d}(pi); q = q_list{d}(pi); m = m_list{d}(pi);
            npxc = size(mov_DNA{p}{q}{m},2);
            npxr = size(mov_DNA{p}{q}{m},1);re
            
            Imax_tmp_DNA = median(cat(1,I_DNA{p}{q}{:}))+1*std(cat(1,I_DNA{p}{q}{:}));
            Imax_tmp_tau = median(cat(1,I_tau{p}{q}{:}))+1*std(cat(1,I_tau{p}{q}{:}));
            Imax_tmp_tau = max(Imax_tmp_tau,20);
        
            ax = subplot(2,4,4*(d-1)+pi); ax.Position(1) = ax.Position(1)+.05 -.02*(pi-1);
            img_DNA = rescale(mov_DNA{p}{q}{m}(:,:,f),0,255, 'InputMin',0,'InputMax',Imax_tmp_DNA);
            img_tau = rescale(mov_tau{p}{q}{m}(:,:,f),0,255, 'InputMin',0,"InputMax",30);
            % img_DNA = img_DNA_avg{p}{q}(:,1:npxc,m);
            % img_tau = img_tau_avg{p}{q}(:,1:npxc,m);
            
            img_DNA_rgb = ind2rgb(round(img_DNA),cmap_DNA);    
            img_tau_rgb = ind2rgb(round(img_tau),cmap_tau);
            img_merge = img_DNA_rgb + img_tau_rgb;
            img_tmp = [img_DNA_rgb; img_tau_rgb; img_merge];
            imshow2(img_tmp,'xdata',[1,npxc]*pxsz,'ydata',[1,3*npxr]*pxsz); hold all;
            xlim([0,111]*pxsz);
            title([legend_list{d}{pi},' tau'],'color','w','fontsize',16);
    
            if pi == 1
                text(-1,3,title_list{d},'color','w','hori','right','fontsize',16);
                if d == 1
                    plot([1,6],[9,9],'w','linewidth',3,'clip','off');
                    text(3.5,9.5,'5 μm','color','w','hori','center','verti','top','fontsize',16);
                end
            end
        end
    end
    text(27,1.2,'DNA','color',hex2rgb(hex_DNA),'fontsize',16,'hori','left');
    text(27,3.4,'Tau','color',hex2rgb(hex_tau),'fontsize',16,'hori','left');
    text(27,5.6,'Merge','color',min(hex2rgb(hex_DNA)+hex2rgb(hex_tau),[1,1,1]),'fontsize',16,'hori','left');
    text(40,3,num2str(f*.1,'t = %.1f s'),'color','w','fontsize',16,'hori','right');
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
end
close(outputVideo);

save('analysis2','sel*');