%% set params
finfo = dir;
pths = {finfo([finfo.isdir]).name}; pths = pths(3:end)'; npth = numel(pths);
tags = {'N100T0','N100T50','N100T500','N100T5000','N25T0','N25T50'};
subpths = {{'raw3'},{'raw1'},{'raw1'},{'raw1'},{'raw1','raw2'},{'raw1','raw2','raw3'}}; nsubpth = cellfun(@(x) numel(x), subpths);

rrange = 1:512; crange = 1:1024; pxsz = 6.5/60;
npxr = numel(rrange); npxc = numel(crange);
threshold = [10*ones(1,4), 20*ones(1,2)];
px_span = 10; npxr_sub = 2*px_span+1; npxc_sub = 150+2*px_span;

varlist = who;
save('analysis',varlist{:});

%% collect data
[fnames,nfile] = deal(cell(size(pths)));
for p = 1:npth
    fnames{p} = cell(nsubpth(p),1); nfile{p} = zeros(nsubpth(p),1);
    [img1_avg,img1_min,img1_max,img2_avg] = deal(cell(nsubpth(p),1));
    for q = 1:nsubpth(p)
        subpth = [pths{p},'\',subpths{p}{q}];
        tifinfo = dir([subpth,'\*.tif']);
        fnames{p}{q} = {tifinfo.name}'; nfile{p}(q) = numel(fnames{p}{q});
        [img1_avg{q},img1_min{q},img1_max{q},img2_avg{q}] = deal(zeros(npxr,npxc,nfile{p}(q)));
        for f = 1:nfile{p}(q)
            disp([p,q,f]);
            mov = loadTifStack16([subpth,'\',fnames{p}{q}{f}]);

            % DNA
            img1_avg{q}(:,:,f) = mean(mov(rrange,crange,:),3);
            img_bkg = medfilt2(img1_avg{q}(:,:,f),[50,50]);
            img1_avg{q}(:,:,f) = img1_avg{q}(:,:,f)-img_bkg;
            img1_min{q}(:,:,f) = double(min(mov(rrange,crange,:),[],3))-img_bkg;
            img1_max{q}(:,:,f) = double(max(mov(rrange,crange,:),[],3))-img_bkg;

            % tau
            img2_avg{q}(:,:,f) = mean(mov(rrange+512,crange,:),3);
            img_bkg = medfilt2(img2_avg{q}(:,:,f),[50,50]);
            img2_avg{q}(:,:,f) = img2_avg{q}(:,:,f)-img_bkg;
        end
    end
    clear mov img_bkg;
    save(['img_',tags{p}],'img*');
end
varlist = union(varlist, {'fnames','nfile'}, 'stable');
save('analysis',varlist{:});

%% calculate channel registration tform
p_list = [1,2,3,6]; q_list = [1,1,1,2];
[movPts_all,fixPts_all,tform] = deal(cell(npth,1));
for p = 1:npth
    load(['img_',tags{p}],'img1_avg','img2_avg');
    [movPts_all{p},fixPts_all{p},tform{p}] = deal(cell(nsubpth(p),1));    
    for q = 1:nsubpth(p)
        if (p == 1 && q == 1) || (p == 2 && q == 1) || (p == 3 && q == 1) || (p == 6 && q == 2)
            [movPts_all{p}{q},fixPts_all{p}{q},tform{p}{q}] = deal([]);
            for f = 1:nfile{p}(q)
                disp([p,q,f]);
                img1_tmp = img1_avg{q}(:,:,f);
                img2_tmp = img2_avg{q}(:,:,f);
    
                img1_tmp = rescale(img1_tmp,"InputMax",median(img1_tmp(:))+4*std(img1_tmp(:)));
                img2_tmp = rescale(img2_tmp,"InputMax",median(img2_tmp(:))+4*std(img2_tmp(:)));
                [movPts,fixPts] = cpselect(img2_tmp,img1_tmp,'Wait',true);
                movPts_all{p}{q} = [movPts_all{p}{q}; movPts];
                fixPts_all{p}{q} = [fixPts_all{p}{q}; fixPts];
            end
            tform{p}{q} = fitgeotform2d(movPts_all{p}{q},fixPts_all{p}{q},"affine");
        end
    end
end
tform{4}{1} = tform{3}{1};
tform{5}{1} = tform{1}{1};
tform{5}{2} = tform{2}{1};
tform{6}{1} = tform{1}{1};
tform{6}{3} = tform{2}{1};

varlist = union(varlist, {'fixPts_all','movPts_all','tform'}, 'stable');
save('analysis',varlist{:});

%% apply channel registration tform
for p = 1:npth
    load(['img_',tags{p}],'img2_avg');
    for q = 1:nsubpth(p)
        for f = 1:nfile{p}(q)
            disp([p,q,f]);
            img2_avg{q}(:,:,f) = imwarp(img2_avg{q}(:,:,f),tform{p}{q},OutputView=imref2d([npxr,npxc]));
        end
    end
    save(['img_',tags{p}],'img2_avg','-append');
end

%% collect points
[nmol,xpos,ypos] = deal(cell(size(pths)));
for p = 1:npth
    load(['img_',tags{p}],'img1*');

    [nmol{p},xpos{p},ypos{p}] = deal(cell(nsubpth(p),1));    
    for q = 1:nsubpth(p)
        subpth = [pths{p},'\',subpths{p}{q}];

        nmol{p}{q} = zeros(nfile{p}(q),1);
        [xpos{p}{q},ypos{p}{q}] = deal(cell(nfile{p}(q),1));
        for f = 1:nfile{p}(q)
            disp([p,q,f]);
            
            maxfig(f); clf; clear h;
            h(1) = subplot(221); imshow3(img1_avg{q}(:,:,f),10); title('avg');
            h(2) = subplot(222); imshow3(img1_min{q}(:,:,f),10); title('min');
            h(3) = subplot(223); imshow3(img1_max{q}(:,:,f),10); title('max');
            h(4) = subplot(224); imshow3(img1_avg{q}(:,:,f),10); hold on;
            linkaxes(h(:),'xy');

            [nmol{p}{q}(f),xpos{p}{q}{f},ypos{p}{q}{f}] = countSM(img1_min{q}(:,:,f),threshold(p),0);
            plot(xpos{p}{q}{f},ypos{p}{q}{f},'k.','markersize',6);
            title([num2str(nmol{p}{q}(f)),' points']);
            
            saveas(gcf,[subpth,'\fullview_',num2str(f),'.jpg']);
        end
    end
end
varlist = union(varlist, {'nmol','xpos','ypos'}, 'stable');
save('analysis',varlist{:});

%% collect dna
sel_m = {{[1:8,10,12:14,19,24,26,28:32,37:38,40:49,51:54,56,59,60,62:63,67:72,74:75,78:82]},...
    {[1:10,12:13,15:16,19,21,23:25]},...
    {[9,12,15,18,20,23,32,36,40:42,45:46,53,56:57,59:61,64:65,67:69,71,73]},...
    {[1,10,12,16,20,25,30,38,57,66]},...
    {[1,6:7,9,13,15:18,21:24,27:29,32:35,38:42,44:47,49:50,52:55,58,61:68,74:75,77:80,84,86,88:89,91:94,96:97,99], [1,3,5,8,10:11,15,17:18,20,21]},...
    {[4,10,52,62:63,65:66,85,93,100,102,104:105,107,112,126,129,144,146,149,154,162,165,168,175,179,181,187,190:192,195,198,200:201,216,219,225,237,240,244,246],[1,5,14,16:17,22,24:25,32,37,39,42,44:45,49,51,55,57,59,62,65,68:69],[4,8,13,17,29,37,43,45,47,54,56:57,62,71,74,80,85,87,93,99]}};

[fk_list,img_DNA_avg,img_DNA_max,img_tau_avg,att_rc,att_dist,nDNA] = deal(cell(npth,1));
for p = 1:npth
    load(['img_',tags{p}],'img*');

    [fk_list{p},img_DNA_avg{p},img_DNA_max{p},img_tau_avg{p},att_rc{p},att_dist{p},nDNA{p}] = deal(cell(nsubpth(p),1));
    for q = 1:nsubpth(p)
        subpth = [pths{p},'\',subpths{p}{q}];
        
        for f = 1:nfile{p}(q)
            disp([p,q,f]);
            [fk_list_tmp,img_DNA_avg_tmp,img_DNA_max_tmp,img_tau_avg_tmp,att_rc_tmp,att_dist_tmp] = deal([]);

            img_tmp1 = img1_avg{q}(:,:,f);
            x_list = xpos{p}{q}{f}; y_list = ypos{p}{q}{f};

            for k1 = 1:nmol{p}{q}(f)
                x1 = x_list(k1); y1 = y_list(k1);
                I1 = mmean(img_tmp1(round(y1)+(-1:1),round(x1)+(-1:1)));
                    
                att_dist_tmp2 = pdist2([x_list,y_list],[x1,y1])*pxsz; % in μm
                k2_list = find(att_dist_tmp2<16 & att_dist_tmp2>3 & x_list>x1); % limit length

                for ki = 1:numel(k2_list)
                    clear test;
                    k2 = k2_list(ki);
                    x2 = x_list(k2); y2 = y_list(k2);
                    I2 = mmean(img_tmp1(round(y2)+(-1:1),round(x2)+(-1:1)));
                    I_anchor = min(I1,I2);

                    pixels = linePixels(x1,y1,x2,y2);
                    r_list = pixels(:,2); c_list = pixels(:,1);
                    idx = sub2ind([npxr,npxc],r_list,c_list);
                    Ivals = img_tmp1(idx);
                    test(1) = all(Ivals > I_anchor*.2);

                    x0 = x1 - (x2-x1)/(att_dist_tmp2(k2)/pxsz)*7;
                    y0 = y1 - (y2-y1)/(att_dist_tmp2(k2)/pxsz)*7;
                    I0 = mmean(img_tmp1(round(y0)+(-1:1),round(x0)+(-1:1)));
                    test(2) = I0 < I1*.5;

                    x3 = x2 + (x2-x1)/(att_dist_tmp2(k2)/pxsz)*7;
                    y3 = y2 + (y2-y1)/(att_dist_tmp2(k2)/pxsz)*7;
                    I3 = mmean(img_tmp1(round(y3)+(-1:1),round(x3)+(-1:1)));
                    test(3) = I3 < I2*.5;

                    I_bkg = min(I0,I3);
                    
                    if all(test)
                        alpha = rad2deg(cart2pol(diff([x1,x2]),diff([y1,y2])));
                        RotMatrix = [cosd(alpha) -sind(alpha); sind(alpha) cosd(alpha)];
                        ori = (([npxr,npxc]+1)/2)';         % Center of the main image
                        for i = 1:3
                            if i == 1
                                img_tmp2 = img1_avg{q}(:,:,f);
                            elseif i == 2
                                img_tmp2 = img1_max{q}(:,:,f);
                            else
                                img_tmp2 = img2_avg{q}(:,:,f);
                            end
                            img_tmp2 = imrotate(img_tmp2,alpha,'bicubic');
                            if i == 1
                                ori_rot = ((size(img_tmp2)+1)/2)';  % Center of the transformed image
                                rc_rot_k1 = RotMatrix*([y1; x1]-ori)+ori_rot;
                                rc_rot_k2 = RotMatrix*([y2; x2]-ori)+ori_rot;
                                r1 = rc_rot_k1(1); c1 = rc_rot_k1(2); c2 = rc_rot_k2(2);
                                rrange_tmp = round(r1)+(-px_span:px_span);
                                crange_tmp = round(c1)-px_span:round(c2)+px_span;
                                
                                r_offset = round(r1)-px_span - 1;
                                c_offset = round(c1)-px_span - 1;
                                r1 = r1-r_offset; c1 = c1-c_offset; c2 = c2-c_offset;
                                r_cnt = round(r1);
                                c_cnt = round(mean([c1,c2]));
                            end
                            img_tmp2 = [img_tmp2(rrange_tmp,crange_tmp), nan(npxr_sub,npxc_sub-numel(crange_tmp))];
                            
                            if i == 1
                                img_DNA_avg_tmp = cat(3,img_DNA_avg_tmp, img_tmp2);
                            elseif i == 2
                                img_DNA_max_tmp = cat(3,img_DNA_max_tmp, img_tmp2);
                            else
                                img_tau_avg_tmp = cat(3,img_tau_avg_tmp, img_tmp2);
                            end
                        end
                        fk_list_tmp = [fk_list_tmp; [f,k1,k2]];
                        att_rc_tmp = [att_rc_tmp; [r1,c1,c2]];
                        att_dist_tmp = [att_dist_tmp; att_dist_tmp2(k2)];
                    end
                end
            end
            nDNA_tmp = size(fk_list_tmp,1);
            
            % clean up grouped
            if nDNA_tmp ~= 0
                rem_m = false(nDNA_tmp,1);
                for m = 1:nDNA_tmp-1
                    if ~rem_m(m)
                        k1 = fk_list_tmp(m,2); k2 = fk_list_tmp(m,3);
                        x1 = x_list(k1); y1 = y_list(k1);
                        x2 = x_list(k2); y2 = y_list(k2);
    
                        idx = find(fk_list_tmp(m+1:end,2)==k1);
                        if ~isempty(idx)
                            k2_test = fk_list_tmp(m+idx,3);
                            x2_test = x_list(k2_test); y2_test = y_list(k2_test);
                            dist_test = pdist2([x2_test,y2_test],[x2,y2]);
                            angle_test = abs(rad2deg(cart2pol(x2_test-x1,y2_test-y1)));
                            idx_rem = find(dist_test < 5 & angle_test < 10);
                            m_list = [m; m+idx(idx_rem)];
                            [~,idx_sel] = max(att_rc_tmp(m_list,1));
                            rem_m(m_list(1:end ~=idx_sel)) = true;
                        end
    
                        idx = find(fk_list_tmp(m+1:end,3)==k2);
                        if ~isempty(idx)
                            k1_test = fk_list_tmp(m+idx,2);
                            x1_test = x_list(k1_test); y1_test = y_list(k1_test);
                            dist_test = pdist2([x1_test,y1_test],[x1,y1]);
                            angle_test = abs(rad2deg(cart2pol(x2-x1_test,y2-y1_test)));
                            idx_rem = find(dist_test < 5 & angle_test < 10);
                            m_list = [m; m+idx(idx_rem)];
                            [~,idx_sel] = max(att_rc_tmp(m_list,1));
                            rem_m(m_list(1:end ~=idx_sel)) = true;
                        end
                    end
                end
                img_DNA_avg{p}{q} = cat(3, img_DNA_avg{p}{q}, img_DNA_avg_tmp(:,:,~rem_m));
                img_DNA_max{p}{q} = cat(3, img_DNA_max{p}{q}, img_DNA_max_tmp(:,:,~rem_m));
                img_tau_avg{p}{q} = cat(3, img_tau_avg{p}{q}, img_tau_avg_tmp(:,:,~rem_m));
                fk_list{p}{q} = [fk_list{p}{q}; fk_list_tmp(~rem_m,:)];
                att_rc{p}{q} = [att_rc{p}{q}; att_rc_tmp(~rem_m,:)];
                att_dist{p}{q} = [att_dist{p}{q}; att_dist_tmp(~rem_m,:)];
            end
        end
        % further discard bad ones
        img_DNA_avg{p}{q} = img_DNA_avg{p}{q}(:,:,sel_m{p}{q});
        img_DNA_max{p}{q} = img_DNA_max{p}{q}(:,:,sel_m{p}{q});
        img_tau_avg{p}{q} = img_tau_avg{p}{q}(:,:,sel_m{p}{q});
        fk_list{p}{q} = fk_list{p}{q}(sel_m{p}{q},:);
        att_rc{p}{q} = att_rc{p}{q}(sel_m{p}{q},:);
        att_dist{p}{q} = att_dist{p}{q}(sel_m{p}{q},:);

        if ~isempty(fk_list{p}{q})
            nDNA{p}{q} = arrayfun(@(f) sum(fk_list{p}{q}(:,1)==f), (1:nfile{p}(q))');
    
            % display results
            for f = nfile{p}(q):-1:1
                img_tmp = img1_avg{q}(:,:,f);
                x_list = xpos{p}{q}{f}; y_list = ypos{p}{q}{f};
                m_list = find(fk_list{p}{q}(:,1) == f);

                maxfig(100+f); clf;
                subplot(121); imshow3(img_tmp,10); title(nDNA{p}{q}(f)); hold on;
                if nDNA{p}{q}(f) ~= 0
                    plot(x_list(fk_list{p}{q}(m_list,2:3))',y_list(fk_list{p}{q}(m_list,2:3))','-');
                    xdat = x_list(fk_list{p}{q}(m_list,2)); ydat = y_list(fk_list{p}{q}(m_list,2));
                    [xdat_uni,ia,ic] = unique(xdat);
                    for xi = 1:numel(xdat_uni)
                        text(xdat_uni(xi),ydat(ia(xi)),num2str(m_list(ic == xi)),'color','k');
                    end
                end

                nrow = max(ceil(nDNA{p}{q}(f)/3),5);
                pid_list = mod(0:nDNA{p}{q}(f)-1,3)+1 + 6*floor((0:nDNA{p}{q}(f)-1)/3) + 3;
                for mi = 1:numel(m_list)
                    m = m_list(mi);
                    subplot(nrow,6,pid_list(mi));
                    img_tmp = img_DNA_max{p}{q}(:,:,m); img_tmp = rescale(img_tmp,"InputMax",nanmedian(img_tmp(:))+2*nanstd(img_tmp(:)));
                    imshow2(img_tmp); xlim([0,npxc_sub]);hold on;
                    plot(att_rc{p}{q}(m,2:3),att_rc{p}{q}(m,1),'b.'); % anchors
                    title(m);
                end
                % pause;
                saveas(100+f,[subpth,'\processed_',num2str(f),'.jpg']);
            end
            close all;
        else
            nDNA{p}{q} = 0;
        end
    end
end
varlist = union(varlist, {'sel_m','fk_list','img_DNA_avg','img_DNA_max','img_tau_avg','att_rc','att_dist','nDNA'}, 'stable');
save('analysis',varlist{:});

%% estimate envelope width
sinefit = fittype('A*sin(x/T*2*pi)','independent','x','problem','T','coeff','A');         
fitOptions = fitoptions('Method','NonlinearLeastSquares','StartPoint', 3);

[mov_DNA,mov_tau,r_maxI,env_fitres,env_gof] = deal(cell(npth,1));
for p = 1:npth
    [mov_DNA{p},mov_tau{p},r_maxI{p},env_fitres{p},env_gof{p}] = deal(cell(nsubpth(p),1));
    for q = 1:nsubpth(p)
        subpth = [pths{p},'\',subpths{p}{q}];
        f = 0;
        [mov_DNA{p}{q},mov_tau{p}{q},r_maxI{p}{q},env_fitres{p}{q},env_gof{p}{q}] = deal(cell(sum(nDNA{p}{q}),1));
        for m = 1:sum(nDNA{p}{q})
            disp([p,q,m]);

            if fk_list{p}{q}(m,1) ~= f
                f = fk_list{p}{q}(m,1);
                mov = loadTifStack16([subpth,'\',fnames{p}{q}{f}]);
                mov1 = double(mov(rrange,crange,:));
                img_bkg = medfilt2(mean(mov1,3),[50,50]);
                mov1 = mov1 - repmat(img_bkg,[1,1,size(mov,3)]);
                
                mov2 = double(mov(rrange+512,crange,:));
                img_bkg = medfilt2(mean(mov2,3),[50,50]);
                mov2 = mov2 - repmat(img_bkg,[1,1,size(mov,3)]);
                mov2 = imwarp(mov2,tform{p}{q},OutputView=imref2d([npxr,npxc]));
            end
    
            k1 = fk_list{p}{q}(m,2); k2 = fk_list{p}{q}(m,3);
            x_list = xpos{p}{q}{f}; y_list = ypos{p}{q}{f};
            x1 = x_list(k1); y1 = y_list(k1);
            x2 = x_list(k2); y2 = y_list(k2);
            alpha = rad2deg(cart2pol(diff([x1,x2]),diff([y1,y2])));
            mov1_rot = imrotate(mov1,alpha,'bicubic');
            mov2_rot = imrotate(mov2,alpha,'bicubic');

            RotMatrix = [cosd(alpha) -sind(alpha); sind(alpha) cosd(alpha)];
            ori = (([npxr,npxc]+1)/2)';         % Center of the main image
            ori_rot = ((size(mov1_rot(:,:,1))+1)/2)';  % Center of the transformed image
            rc_rot_k1 = RotMatrix*([y1; x1]-ori)+ori_rot;
            rc_rot_k2 = RotMatrix*([y2; x2]-ori)+ori_rot;

            r1 = rc_rot_k1(1); c1 = rc_rot_k1(2); c2 = rc_rot_k2(2);
            rrange_tmp = round(r1)+(-px_span:px_span);
            crange_tmp = round(c1)-px_span:round(c2)+px_span;

            mov_DNA{p}{q}{m} = uint16(mov1_rot(rrange_tmp,crange_tmp,:));
            mov_tau{p}{q}{m} = uint16(mov2_rot(rrange_tmp,crange_tmp,:));
            
            mov_DNA_tmp = mov1_rot(rrange_tmp,crange_tmp,:);

            filterSize = 5; kernel = ones(filterSize) / filterSize^2;
            tmp = arrayfun(@(t) conv2(mov_DNA_tmp(:,:,t), kernel, 'same'), 1:size(mov_DNA_tmp,3),'unif',0);
            mov_DNA_tmp = cat(3,tmp{:});
            [~,r_maxI{p}{q}{m}] = max(mov_DNA_tmp,[],1);
            idx1 = squeeze(max(r_maxI{p}{q}{m},[],3));
            idx2 = squeeze(min(r_maxI{p}{q}{m},[],3));
            
            c1 = att_rc{p}{q}(m,2); c2 = att_rc{p}{q}(m,3); offset = 5;
            crange_tmp = round(c1)+offset:round(c2)-offset;
            xdat_fit = (0:numel(crange_tmp)-1)+offset;
            ydat1 = idx1(crange_tmp)-(px_span+1);
            ydat2 = -(idx2(crange_tmp)-(px_span+1));
            ydat_fit = mean([ydat1;ydat2]);
            [env_fitres{p}{q}{m},env_gof{p}{q}{m}] = fit(xdat_fit',ydat_fit',sinefit,'problem',(round(c2)-round(c1))*2,fitOptions);
        end
    end
end
clear mov1* mov2* *tmp;
save('movies','mov_DNA','mov_tau');
varlist = union(varlist, {'sinefit','fitOptions','r_maxI','env_fitres','env_gof'}, 'stable');
save('analysis',varlist{:});

%% extract intensity and correlation
rrange_tmp = (px_span+1)+(-1:1);
[I_DNA,I_tau,RMSd,pearson] = deal(cell(npth,1));
for p = 1:npth
    [I_DNA{p},I_tau{p},RMSd{p},pearson{p}] = deal(cell(nsubpth(p),1));
    for q = 1:nsubpth(p)
        nDNA_tmp = sum(nDNA{p}{q});
        [I_DNA{p}{q},I_tau{p}{q}] = deal(cell(nDNA_tmp,1));
        [RMSd{p}{q},pearson{p}{q}] = deal(nan(nDNA_tmp,2));
        for m = 1:nDNA_tmp
            c1 = att_rc{p}{q}(m,2); c2 = att_rc{p}{q}(m,3); offset = 0;
            crange_tmp = round(c1)+offset:round(c2)-offset;
            I_DNA_tmp = mean(img_DNA_avg{p}{q}(rrange_tmp,crange_tmp,m),1)';
            I_tau_tmp = mean(img_tau_avg{p}{q}(rrange_tmp,crange_tmp,m),1)';
            if max(I_DNA_tmp)>1600
                disp([p,q,m]);
            else
                I_DNA_tmp(I_DNA_tmp<0) = 0;
                
                I_tau_tmp = I_tau_tmp-I_DNA_tmp*crosstalk{p}(q); % crosstalk correction
                I_tau_tmp(I_tau_tmp<0) = 0;
                
                I_DNA{p}{q}{m} = I_DNA_tmp;
                I_tau{p}{q}{m} = I_tau_tmp;
                RMSd{p}{q}(m,1) = rms(diff(I_DNA_tmp))/mean(I_DNA_tmp);
                RMSd{p}{q}(m,2) = rms(diff(I_tau_tmp))/mean(I_tau_tmp);
                pearson{p}{q}(m,1) = corr(I_DNA_tmp,I_tau_tmp);
                pearson{p}{q}(m,2) = corr(I_DNA_tmp,flipud(I_tau_tmp));
                % pearson{p}{q}(m,1) = corr(diff(I_DNA_tmp),diff(I_tau_tmp));
                % pearson{p}{q}(m,2) = corr(diff(I_DNA_tmp),flipud(diff(I_tau_tmp)));
            end
        end
    end
end
varlist = union(varlist, {'crosstalk','I_DNA','I_tau','RMSd','pearson'}, 'stable');
save('analysis',varlist{:});

%% display results
% load('analysis');
crosstalk = {.0416, .0422, .0536, .0536, [.0445,.0378], [.0445,.0435,.0378]};
for p = 1:npth
    for q = 1:nsubpth(p)
        Imax_tmp_DNA = median(cat(1,I_DNA{p}{q}{:}))+1*std(cat(1,I_DNA{p}{q}{:}));
        Imax_tmp_tau = median(cat(1,I_tau{p}{q}{:}))+1*std(cat(1,I_tau{p}{q}{:}));
        Imax_tmp_tau = max(Imax_tmp_tau,10);

        maxfig(100+10*p+q); clf; clear ax;
        nrow = max(ceil(sum(nDNA{p}{q})/6),5);
        for m = 1:sum(nDNA{p}{q})
            rid = floor((m-1)/6); cid = mod(m-1,6);
            ax(m) = axes('Position',[(1/6)*cid+.02,(nrow-rid-1)*(1/nrow),1/6-.02,(1/nrow)*.8]);
            
            img_tmp1 = img_DNA_max{p}{q}(:,:,m);
            img_tmp2 = img_DNA_avg{p}{q}(:,:,m);
            img_tmp3 = img_tau_avg{p}{q}(:,:,m); % avg img of tau
            img_tmp3 = img_tmp3 - img_tmp2*crosstalk{p}(q); % crosstalk correction

            img_tmp1 = rescale(img_tmp1,'InputMin',0,"InputMax",Imax_tmp_DNA);
            img_tmp2 = rescale(img_tmp2,'InputMin',0,"InputMax",Imax_tmp_DNA);
            img_tmp3 = rescale(img_tmp3,'InputMin',0,"InputMax",Imax_tmp_tau);
            
            img_tmp = [img_tmp1; img_tmp2; img_tmp3];
            imshow2(img_tmp); axis manual; xlim([0,npxc_sub]); hold on;

            c1 = att_rc{p}{q}(m,2); c2 = att_rc{p}{q}(m,3); offset = 5;
            crange_tmp = round(c1)+offset:round(c2)-offset;
            xdat_fit = (0:numel(crange_tmp)-1)+offset;
            idx1 = squeeze(max(r_maxI{p}{q}{m},[],3));
            idx2 = squeeze(min(r_maxI{p}{q}{m},[],3));            
            plot(crange_tmp,idx1(crange_tmp),'r.','markersize',5);
            plot(crange_tmp,idx2(crange_tmp),'g.','markersize',5);
            if ~isnan(env_fitres{p}{q}{m}.A)
                plot(crange_tmp,env_fitres{p}{q}{m}(xdat_fit)+px_span+1,'y');
                plot(crange_tmp,-env_fitres{p}{q}{m}(xdat_fit)+px_span+1,'y');
            end
            text(0,npxr_sub/2,num2str(m),'hori','right');
            title(num2str([env_fitres{p}{q}{m}.A*2*pxsz,env_gof{p}{q}{m}.rmse],'env = %.1f μm, rmse = %.2f'));
        end
        % beep; pause;
        saveas(100+10*p+q,['collected_',tags{p},'_',subpths{p}{q},'.jpg']);
    end
end
% something wrong (tau image clipped) in 6,1,8

%% export movies
for p = 1:npth
    for q = 1:nsubpth(p)
        subpth = [pths{p},'\',subpths{p}{q}];
        Imax_tmp_DNA = median(cat(1,I_DNA{p}{q}{:}))+1*std(cat(1,I_DNA{p}{q}{:}));
        Imax_tmp_tau = median(cat(1,I_tau{p}{q}{:}))+1*std(cat(1,I_tau{p}{q}{:}));
        Imax_tmp_tau = max(Imax_tmp_tau,10);
        
        for m = 1:sum(nDNA{p}{q})
            disp([p,q,m]);
            outputVideo = VideoWriter([subpth,'\DNA_',num2str(m),'.avi'], 'Grayscale AVI');
            outputVideo.FrameRate = 10;
            open(outputVideo);
            for i = 1:size(mov_DNA{p}{q}{m}, 3)
                tmp = uint8(rescale(mov_DNA{p}{q}{m}(:,:,i),'InputMin',0,'InputMax',Imax_tmp_DNA)*255);
                writeVideo(outputVideo, tmp);
            end
            close(outputVideo);

            outputVideo = VideoWriter([subpth,'\tau_',num2str(m),'.avi'], 'Grayscale AVI');
            outputVideo.FrameRate = 10;
            open(outputVideo);
            for i = 1:size(mov_tau{p}{q}{m}, 3)
                tmp = uint8(rescale(mov_tau{p}{q}{m}(:,:,i),'InputMin',0,'InputMax',Imax_tmp_tau)*255);
                writeVideo(outputVideo, tmp);
            end
            close(outputVideo);
            % pause;
        end
    end
end

%% bleaching check
[I_blc] = cell(npth,1);
maxfig(101); clf;
pid_list = {1, 2, 3, 3, [1,2], [1,4,2]};
title_list = {'240106','230915','231110','230921'};
expfit = fittype('I1*exp(-(f*.1/tau))+I0','indep','f','coeff',{'tau','I1','I0'});
for p = 1:npth
    I_blc{p} = zeros(50,nsubpth(p),2);
    for q = 1:nsubpth(p)
        I_tmp = cellfun(@(x) squeeze(mmean(x,[1,2])), mov_DNA{p}{q}, 'unif',0);
        I_blc{p}(:,q,1) = mean(cat(2,I_tmp{:}),2);
        subplot(2,4,pid_list{p}(q));
        h = plot(I_blc{p}(:,q,1)); hold all;
        title([title_list{pid_list{p}(q)},', DNA']);
        xdat_fit = (0:(size(I_blc{p},1)-1))';
        ydat_fit = I_blc{p}(:,q,1);
        fitres = fit(xdat_fit,ydat_fit,expfit,'start',[2,ydat_fit(1)-ydat_fit(end),ydat_fit(end)]);
        plot(xdat_fit,fitres(xdat_fit),'k:');
        disp(fitres.I0);

        I_tmp = cellfun(@(x) squeeze(mmean(x,[1,2])), mov_tau{p}{q}, 'unif',0);
        I_blc{p}(:,q,2) = mean(cat(2,I_tmp{:}),2);
        subplot(2,4,pid_list{p}(q)+4);
        plot(I_blc{p}(:,q,2)); hold all;
        title([title_list{pid_list{p}(q)},', tau']);
        ydat_fit = I_blc{p}(:,q,2);
        fitres = fit(xdat_fit,ydat_fit,expfit,'start',[2,ydat_fit(1)-ydat_fit(end),ydat_fit(end)]);
        plot(xdat_fit,fitres(xdat_fit),'k:');
        disp(fitres.I0);
    end
end

%% kymographs
% load('movies');
% 6,1,41: interference from nearby DNA
for p = 1:npth
    for q = 1:nsubpth(p)
        maxfig(200+10*p+q); clf; clear ax;
        nrow = max(ceil(sum(nDNA{p}{q})/6),5);
        Imax_tmp_DNA = median(cat(1,I_DNA{p}{q}{:}))+1*std(cat(1,I_DNA{p}{q}{:}));
        Imax_tmp_tau = median(cat(1,I_tau{p}{q}{:}))+1*std(cat(1,I_tau{p}{q}{:}));
        Imax_tmp_tau = max(Imax_tmp_tau,10);
        for m = 1:sum(nDNA{p}{q})
            rid = floor((m-1)/6); cid = mod(m-1,6);
            
            mov_tmp1 = mov_DNA{p}{q}{m};
            mov_tmp2 = mov_tau{p}{q}{m};
            
            rrange_tmp = (px_span+1)+(-5:5);
            c1 = att_rc{p}{q}(m,2); c2 = att_rc{p}{q}(m,3); offset = -7;
            crange_tmp = round(c1)+offset:round(c2)-offset;
            img_tmp1 = squeeze(mean(mov_tmp1(rrange_tmp,crange_tmp,:),1))';
            img_tmp1 = rescale(img_tmp1,'InputMin',0,'InputMax',Imax_tmp_DNA);
        
            img_tmp2 = squeeze(mean(mov_tmp2(rrange_tmp,crange_tmp,:),1))';
            img_tmp2 = rescale(img_tmp2,'InputMin',0,'InputMax',Imax_tmp_tau);
            img_kym = cat(3,img_tmp2,img_tmp1,zeros(size(img_tmp1)));
            
            ax(m) = axes('Position',[(1/6)*cid+.02,(nrow-rid-1)*(1/nrow)+.01,(1/6)*(numel(crange_tmp)/140),(1/nrow)*.8]);
            imshow2(img_kym);
            title(m);
        end
        saveas(200+10*p+q,['kym_',tags{p},'_',subpths{p}{q},'.jpg']);
    end
end
