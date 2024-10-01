b_list = [2,2,1; 1,1,1; 2,2,2; 1,1,1; 1,1,1];
f_list = zeros(npth,3,3,2); % nfile, repeat, direction
for n = 1:nfile(p)
    for p = 1:npth
        if p == 1
            maxfig(n); clf;
        else
            figure(n);
        end

        subplot(5,1,p);
        b = b_list(p,n); return_dat;

        f0 = 1;
        for i = 1:3
            f_list(p,n,i,1) = find(Fdat(f0:end)<2,1,'first') +(f0-1);
            f0 = f_list(p,n,i,1) + 100;
            f_list(p,n,i,2) = find(Fdat(f0:end)>2,1,'first') +(f0-1);
            f0 = f_list(p,n,i,2) + 100;
        end
            
        % frange = 1:nframe{p}(n);
        frange = f_list(p,n,1,1)-500:f_list(p,n,3,2)+500;
        yyaxis left;
        plot(tdat(frange),Fdat(frange)); hold all;
        vline(reshape(tdat(f_list(p,n,:,1)),[],1),'b:');
        vline(reshape(tdat(f_list(p,n,:,2)),[],1),'r:');
        axis tight;
        yyaxis right;
        plot(tdat(frange),dzdat(frange));
        ylim([-3e3,300]);
        title([p,n]);
        % pause;
    end
end 

%%
maxfig(11); clf;
for p = 1:npth
    for n = 1:3
        b = b_list(p,n); return_dat;
        for j = 1:2
            subplot(5,6,6*(p-1)+3*(j-1)+n);
            for i = 1:3
                frange = f_list(p,n,i,j)+(-200:3000);
                tdat_tmp = tdat(frange)-tdat(frange(1));
                dzdat_tmp = dzdat(frange);
                dzdat_tmp = medfilt1(dzdat_tmp,11);
                plot(tdat_tmp,dzdat_tmp); hold all;
            end
            xlim([0,30]);
            if j == 1
                ylim([-3e3,300]);
            else
                ylim([-500,0]);
            end
        end
    end
end

%%
maxfig(12); clf;
[tdat_all,dzdat_all,ddzdat_all] = deal(cell(2,3));
offset = 2;
for j = 1:2
    for n = 1:3
        h(1) = subplot(4,3,6*(j-1)+n);
        for p = 1:npth
            b = b_list(p,n); return_dat;
            % dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n}(b),ori{p}{n}(b));
            dzdat = dzdat-median(dzdat(abs(Fdat-5)<.1));
            for i = 1:3
                frange = f_list(p,n,i,j)+(-200:3000);
                tdat_tmp = tdat(frange);
                dzdat_tmp = dzdat(frange);
                dzdat_tmp = medfilt1(dzdat_tmp,21);
                plot(tdat_tmp-tdat_tmp(1),dzdat_tmp,'linewidth',1); hold all;
                ddzdat_tmp = dzdat_tmp((offset+1):end)-dzdat_tmp(1:end-offset);
                ddzdat_tmp = smooth(ddzdat_tmp,7);


                tdat_all{j,n} = [tdat_all{j,n}, tdat_tmp];
                dzdat_all{j,n} = [dzdat_all{j,n}, dzdat_tmp];
                ddzdat_all{j,n} = [ddzdat_all{j,n}, ddzdat_tmp];
            end
        end
        xlim([0,numel(frange)*.01]);
        if j == 1
            ylim([-3e3,300]);
        else
            ylim([-500,0]);
        end
        
        h(2) = subplot(4,3,6*(j-1)+n+3);
        plot(tdat_tmp(1:end-offset)-tdat_tmp(1),ddzdat_all{j,n});
        xlim([0,numel(frange)*.01]);
        linkaxes(h,'x');
    end
end

%%
stepsize = cell(3,1);
j = 2;
for n = 1:3
    maxfig(20+n); clf;
    for p = 1:15
        subplot(3,5,p);
        plot(dzdat_all{j,n}(:,p));
        ylim([-500,0]);
        [pks,locs] = findpeaks(ddzdat_all{j,n}(:,p),'MinPeakProminence',10);
        % locs(locs<206) = [];
        locs(locs<221) = [];
        if any(locs)
            vline(locs);
        end

        for i = 1:numel(locs)
            % stepsize{n} = [stepsize{n}; median(dzdat_all{j,n}(locs(i)+(3:6),p)) - median(dzdat_all{j,n}(locs(i)+(-4:-1),p))];
            stepsize{n} = [stepsize{n}; mean(dzdat_all{j,n}(locs(i)+(3:5),p)) - mean(dzdat_all{j,n}(locs(i)+(-3:-1),p))];
        end
    end
    
    figure(31);
    subplot(2,3,3*(j-1)+n);
    % histogram(stepsize{n});
    histogram(stepsize{n},0:25:250);
    vline(median(stepsize{n}));
    title(numel(stepsize{n}));

    if n == 1
        locs
    end
end

save('analysis2_FC','*_list','*dat_all','stepsize');