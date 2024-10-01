function co_pos=colocal_1d(fname,img_1,img_2,feature_thres,max_off, disp_range)
    %img_1 for upper, img_2 for below
    if numel(size(img_1))==3 nframe=size(img_1,3);
    else nframe=1; end

    for imageNumber=1:nframe
        if imageNumber==1 draw=1;
        else draw=0; end
        [nmol1,xpos1,ypos1]=countSM(img_1(:,:,imageNumber),feature_thres(1),draw);

        if nmol1==0 
            co_pos{imageNumber}=[]; 
            if imageNumber==1 co_ratio(imageNumber,1)=0; end
            continue;
        end

        temp_pos=[];
        for n=1:nmol1
            roi=movmean(img_2(round(ypos1(n)),round(xpos1(n)-5:xpos1(n)+5),imageNumber),[3 3]);
            if abs(find(roi==max(roi))-6)<max_off
            temp_pos=[temp_pos; xpos1(n) ypos1(n)]; 
            end
        end
        co_pos{imageNumber}=temp_pos;
    end

   %% comp image
    figure;   
   
    img_comp=double(cat(3, img_2(:,:,1), img_1(:,:,1), zeros(size(img_2(:,:,1))))); 
    img_comp(:,:,1)=img_comp(:,:,1)/disp_range(2); img_comp(:,:,2)=img_comp(:,:,2)/disp_range(1);
    imshow2(img_comp); hold on;
    title(fname);
    if ~isempty(co_pos{1})
        plot(co_pos{1}(:,1),co_pos{1}(:,2),'yo');
    end
    hold off
    saveas(gcf,fname)
end