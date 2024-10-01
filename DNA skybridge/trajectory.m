function out=trajectory(sz, tr)
%function for calculation trajectories of particles
%input:
%   tr: result of the function "track"
%       (x y t id) matrix, which is sorted according to the time
%output:
%   out: all informations of lines
%        (xi; yi; xf; yf; t; id)

ppre=[];    %x y t id
img=zeros([sz(1) sz(2)]);
    for imageNumber = 1:sz(3)
        [r,~]=find(tr(:,3)==imageNumber);  %t초에 해당하는 frame에 존재하는 particle rows
        line=[];
        for p=1:length(r)   %particle
            if isempty(ppre) || isempty(find(ppre(:,4)==tr(r(p),4)))    %처음 발견한 particle
                ppre=vertcat(ppre, tr(r(p),:));
            else 
                [rp, ~]=find(ppre(:,4)==tr(r(p),4));   %다시 발견한 particle
                [x,y]=bresenham(ppre(rp, 1),ppre(rp, 2), tr(r(p),1), tr(r(p),2));
                img(sub2ind(size(img),y,x))=1;
                ppre(rp,:)=tr(r(p),:);
            end
        end
        out(:,:,imageNumber)=img;
%         comp_img=cat(3,a(:,:,imageNumber),img*255,zeros(size(a(:,:,imageNumber))));
%         writeTIFrgb(comp_img,fname);
    end
end

