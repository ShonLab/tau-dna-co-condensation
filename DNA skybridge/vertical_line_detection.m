function out=vertical_line_detection(img, barwidth,period)
    bkg=medfilt2(img,round(size(img)/10));
    img1=medfilt2(img,[size(img,1),barwidth/2]);
    SE = strel("line",size(img,1),90);
    filtered=imdilate(img1,SE);

%     mask=imbinarize(filtered,'adaptive');
%     mask=~imdilate(mask,strel("line",barwidth,0));

    mask=~imbinarize(filtered,'adaptive');

    img2=img.*uint16(mask);
    T=adaptthresh(img2,0.5,"Statistic","gaussian");
    img2=imbinarize(img2,T);
    img2=imdilate(img2,strel("line",3,0));
    

%     img2=imbinarize(medfilt2(img.*uint16(mask)+bkg.*uint16(~mask),[1,3]),'adaptive');
%     img2=img2.*mask;
    mask2=zeros(size(img));

    height=1;

    for row=1:size(img,1)
        run_length = img2(row,:);
        count=0;
        for idx=1:size(img,2)
            if run_length(idx)==1
                count=count+1; 
                if idx==size(img,2)
                    if count>=period
                        mask2(max(1,row-height):min(row+height,size(img,1)),idx-count:idx-1)=1;
                    end
                end
            else
                if count>=period
                    mask2(max(1,row-height):min(row+2,size(img,1)),idx-count:idx-height)=1;
                end
                count=0;
            end
        end

    end

    mask=mask & mask2;
    out=mask;
end