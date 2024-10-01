function [img_1,img_2]=dualviewer_merger(movie,tform,fname)
    nframe=size(movie,numel(size(movie)));
    delete(fname)

    up=movie(1:end/2,:,:);
    down=movie(1+end/2:end,:,:);     

    for imageNumber=1:nframe
        img_1(:,:,imageNumber)=up(:,:,imageNumber);
        img_2(:,:,imageNumber)=down(:,:,imageNumber);

        img_1(:,:,imageNumber) = imwarp(img_1(:,:,imageNumber), tform, 'OutputView', imref2d(size(img_2)));
%         comp_img=double(cat(3, img_2, img_1, zeros(size(img_2))));
%         comp_img=cat(3, img_2(:,:,imageNumber), img_1(:,:,imageNumber), zeros(size(img_2(:,:,imageNumber))));
%         writeTIFrgb(comp_img,fname);        
    end
end
