function out=kymo(roi)
    for frame=1:size(roi,3)
        out(frame,:)=roi(:,:,frame);
    end
end