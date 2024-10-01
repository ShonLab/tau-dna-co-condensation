%% adjustable
%experimental condition
% framerate=0.03;  %second
% scale=216*10^(-9);  %meter
% sz=7; %diameter of particle, unit:pixel
% 
% ball_radius=20;
% threshold=15; dispres=0; 
% 
% len=10; %얼마나 오랜 프레임동안 유지되어야 계산할지
% maxdisp=7; %the maximum distance particles can move after one frame
% param.mem=2; param.dim=2; param.good=1; param.quiet=1; %parameters used in track.m file

% path={'raw1','raw2'};

function out=particle_tracking(img_filtered,threshold,maxdisp,len) 
    out=[]; last_id=0; pre=[];
    for imageNumber=1:size(img_filtered,3)        
        [nmol,xpos,ypos] = countSM(img_filtered(:,:,imageNumber),threshold,0);

        if ~nmol pre=[]; continue;  end
        if isempty(pre)

            pre=[xpos, ypos, repmat(imageNumber,numel(xpos),1), (last_id+1:last_id+numel(xpos))'];
            last_id=last_id+numel(xpos);

        else

            temp=[xpos, ypos]; cur=[];
            for idx=1:numel(xpos)
                cost=temp(idx,:)-pre(:,1:2); cost=sqrt(cost(:,1).^2+cost(:,2).^2);
                if min(cost)<=maxdisp
                    id=find(cost==min(cost),1);
                    cur=[cur; temp(idx,1) temp(idx,2) imageNumber pre(id,4)];
                elseif size(img_filtered,3)-imageNumber>len
                    last_id=last_id+1;
                    cur=[cur; temp(idx,1) temp(idx,2) imageNumber last_id];
                end
            end

            out=[out;cur];

            if isempty(cur) fin=pre(:,4);
            else fin=setdiff(pre(:,4),cur(:,4)); end

            if ~isempty(fin) 
                for f=1:numel(fin)
                if numel(find(out(:,4)==fin(f)))<len
                    out(find(out(:,4)==fin(f)),:)=[];
                end
                end
            end

            pre=cur;
        end
    end
    
    if ~isempty(cur)
        for f=1:numel(cur(:,4))
            if numel(find(out(:,4)==cur(f,4)))<len
                out(find(out(:,4)==cur(f,4)),:)=[];
            end
        end
    end


end




