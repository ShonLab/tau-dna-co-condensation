function out=CalD_1D(dim, tr, framerate, scale,fname)
% calculate the diffusion coefficient using the formula:
% msd=2nDt-->D=msd/2nt
%input:
%   dim: the dimension of the area where particles is moving in
%   lineM: line matrix, the result of the function "trajectory"
%          (x; y; t; id)
%   framerate: in unit 'sec'
%   scale: the pixel size, in unit 'meter'
%
%output:
%   out: the list of the estimated diffusion coefficients of the each particle
%        (D id)
%        in unit m^2/s


    tr=sortrows(tr,4);    %id에 따라 분류
    particles=unique(tr(:,end));  %particle list
    msd={}; %time, msd
    

    for p=1:size(particles,1)
        [r,~]=find(tr(:,4)==particles(p));
%         if size(r,1)<len; continue; end
        for frame=1:size(r,1)-1
            imageNumber=tr(r(frame),3);
            msd{p}{imageNumber,1}=imageNumber*framerate; 
            msd{p}{imageNumber,2}=((tr(r(frame),1)-tr(r(frame)+1,1))^2)*scale^2;
        end            
    end

    %% check if it is normal diffusion
%     figure();
%     D=[];
%     for p=1:size(msd,2)
%         if isempty(msd{p}) continue; end
%         x=cell2mat(msd{p}(:,1)); tmp_y=cell2mat(msd{p}(:,2));
%         y=[];
%         for frames=1:size(tmp_y,1)
%             y(frames)=sum(tmp_y(1:frames));
%         end
%         plot(x,y*10^12); hold on;
%         D=vertcat(D,[(y(end)-y(1))/2/dim/(x(end)-x(1)), particles(p)]);
%     end
%     hold off;
% 
%     ylabel('msd (μm^2)')
%     xlabel('time (s)')
%     title('msd cumulation')
%     saveas(gcf,strcat(fname,'_msd.png'))
%     
%     out=D;

%% dt vs msd plot
    figure();
    D=[]; 
    for p=1:size(msd,2)
        if isempty(msd{p}) continue; end
        tmp_y=cell2mat(msd{p}(:,2));
        y=step_sum(tmp_y);
        y=log10(y); dt=log10((1:size(y,1))*framerate);


        plot(dt,y+12); hold on;
        coeffs=polyfit(dt,y,1);

        D=vertcat(D,[10^coeffs(2)/2/dim, coeffs(1), particles(p)]);
    end
    hold off;

    ylabel('log msd (μm^2)')
    xlabel('log Δt (s)')
    title(fname)
    saveas(gcf,fname)
    
    out=D;

end

function out=step_sum(A)
    n=size(A,1);
    out=[];
    for step=1:n
        tmp=[];
        for idx=1:size(A,1)-step+1
            tmp=vertcat(tmp,sum(A(idx:idx+step-1)));
        end
        out=vertcat(out,mean(tmp));
    end
end
