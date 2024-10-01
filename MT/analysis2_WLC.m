b_list = [2,2,1; 1,1,1; 2,2,2; 1,1,1; 1,1,1];
f_list = zeros(npth,3,3,2,2); % nfile, repeat, direction, start/end
maxfig(1); clf;
for p = 1:npth
    for n = 1:nfile(p)
        subplot(3,5,5*(n-1)+p);
        b = b_list(p,n); return_dat;

        for i = 1:3
            f_tmp = (i-1)*280e2;
            frange_tmp = max(f_tmp-50e2,1):(f_tmp+50e2);
            f_list(p,n,i,1,1) = find(Fdat(frange_tmp)<.05,1,'last') +(frange_tmp(1)-1);

            f_tmp = i*280e2;
            frange_tmp = (f_tmp-50e2):min(f_tmp+50e2,nframe{p}(n));
            f_list(p,n,i,2,2) = find(Fdat(frange_tmp)<.05,1,'first') +(frange_tmp(1)-1);

            frange_tmp = f_list(p,n,i,1,1):f_list(p,n,i,2,2);
            Fdat_tmp = Fdat(frange_tmp);
            f_middle = round(mean(find(Fdat_tmp == max(Fdat_tmp)))) +(frange_tmp(1)-1);            
            f_list(p,n,i,1,2) = f_middle;
            f_list(p,n,i,2,1) = f_middle;
        end
            
        % frange = 1:nframe{p}(n);
        frange = f_list(p,n,1,1,1):f_list(p,n,3,2,2);
        yyaxis left;
        plot(tdat(frange),Fdat(frange)); hold all;
        vline(reshape(tdat(f_list(p,n,:,:,1)),[],1),'b:');
        vline(reshape(tdat(f_list(p,n,:,:,2)),[],1),'r:');
        yyaxis right;
        plot(tdat(frange),dzdat(frange));
        title([p,n]);
        % pause;
    end
end 

%%
maxfig(2); clf;
colororder = colororder;
Fdat_model = logspace(-1,1,100)';
zdat_model = WLC_inv(Fdat_model,10e3*.338,45,300,1)/1e3;
Fdat_ref = .1:.1:10;
[F_unzip,F_zip] = deal(nan(npth,3,3)); % tau conc, repeat
for p = 1:npth
    for n = 1:nfile(p)
        subplot(3,5,5*(n-1)+p);
        b = b_list(p,n); return_dat;
        dzdat2 = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n}(b),ori{p}{n}(b));
        dzdat2 = dzdat2 - mean(dzdat2(abs(Fdat-10)<.1)) + WLC_inv(10,10e3*.338,45,300,1);
        dzdat2 = medfilt1(dzdat2,11)/1e3;

        if n == 1 % construct reference FEC
            frange = f_list(p,n,1,1,1):f_list(p,n,3,2,2);
            dzdat_ref = zeros(size(Fdat_ref));
            for Fi = 1:numel(Fdat_ref)
                grange = find(abs(Fdat(frange)-Fdat_ref(Fi))<.05);
                dzdat_ref(Fi) = mean(dzdat2(frange(grange)));
            end
            dzdat_ref = dzdat_ref*.95;
        end
        
        for i = 1:3
            % stretching
            frange = f_list(p,n,i,1,1):f_list(p,n,i,1,2);
            semilogy(dzdat2(frange),Fdat(frange),'color',colororder(2,:),'linewidth',1); hold all;
            [~,idx] = closest(Fdat_model,Fdat(frange));
            % fi_cross = find(dzdat2(frange) < zdat_model(idx)*.9, 1, 'last');
            fi_cross = find(dzdat2(frange) < interp1(Fdat_ref,dzdat_ref,Fdat(frange)), 1, 'last');
            if ~isempty(fi_cross)
                hline(Fdat(frange(fi_cross)),'r--');
                F_unzip(p,n,i) = Fdat(frange(fi_cross));
            end

            % condensation
            frange = f_list(p,n,i,2,1):f_list(p,n,i,2,2);
            semilogy(dzdat2(frange),Fdat(frange),'color',colororder(1,:),'linewidth',1);
            [~,idx] = closest(Fdat_model,Fdat(frange));
            % fi_cross = find(dzdat2(frange) < zdat_model(idx)*.9, 1, 'first');
            fi_cross = find(dzdat2(frange) < interp1(Fdat_ref,dzdat_ref,Fdat(frange)), 1, 'first');
            if ~isempty(fi_cross)
                hline(Fdat(frange(fi_cross)),'b--');
                F_zip(p,n,i) = Fdat(frange(fi_cross));
            end
        end
        semilogy(zdat_model,Fdat_model,'k--');
        semilogy(dzdat_ref,Fdat_ref,'g--');
        xlim([-.5,3.5]); ylim([.1,10]);
        title([p,n]);
    end
end

%%
maxfig(3); clf;
title_list = {'No tau','500 nM tau','5 μM tau'};
for n = 1:3
    subplot(2,4,n);
    ydat_all = squeeze(F_unzip(:,n,:));
    ydat_avg = mean(ydat_all,2);
    ydat_err = std(ydat_all,[],2);
    bar(ydat_avg,'facecolor','r'); hold all;
    errorbar(ydat_avg,ydat_err,'k','linestyle','none');
    ylim([0,4]);
    title(title_list{n});

    subplot(2,4,n+4);
    ydat_all = squeeze(F_zip(:,n,:));
    ydat_avg = mean(ydat_all,2);
    ydat_err = std(ydat_all,[],2);
    bar(ydat_avg,'facecolor','b'); hold all;
    errorbar(ydat_avg,ydat_err,'k','linestyle','none');
    ylim([0,4]);
    title(title_list{n});
end

subplot(244);
ydat_all = reshape(permute(F_unzip(:,2:3,:),[1,3,2]),15,[]);
ydat_avg = nanmean(ydat_all,1);
ydat_err = nanstd(ydat_all,[],1);
bar(ydat_avg,'facecolor','r'); hold all;
errorbar(ydat_avg,ydat_err,'k','linestyle','none');

subplot(248);
ydat_all = reshape(permute(F_zip(:,2:3,:),[1,3,2]),15,[]);
ydat_avg = nanmean(ydat_all,1);
ydat_err = nanstd(ydat_all,[],1);
bar(ydat_avg,'facecolor','b'); hold all;
errorbar(ydat_avg,ydat_err,'k','linestyle','none');

% save('analysis2_WLC','f_list','F_unzip','F_zip');