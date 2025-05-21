function dat = loadTifStack16(fname,varargin)
% loadTifStack16(fname): Read all frames
% loadTifStack16(fname,nframe): Read nframe frames from the beginning
% loadTifStack16(fname,nframe,f0): Read nframe frames from the f0-th frame

try
    TifLink = Tiff(fname, 'r');
    exitflag = 0;
catch
    disp([fname,' does not exist!']);
    exitflag = 1;
end

if ~exitflag
    % image dimension
    npxr = TifLink.getTag('ImageLength');
    npxc = TifLink.getTag('ImageWidth');
    % frame index to start reading
    if nargin == 3
        f0 = varargin{2};
    else
        f0 = 1;
    end
    % number of frames to read
    if nargin >= 2
        nframe = varargin{1};
    else
        info = imfinfo(fname);
        nframe = numel(info);
    end
    
    % Read frames
    dat = zeros(npxr,npxc,nframe,'uint16');
    for j = 1:nframe
       try
           TifLink.setDirectory(f0-1+j);
       catch
            dat = dat(:,:,1:j-1);
            disp(['Warning: Only ',num2str(j-1),' frames were loaded.']);           
           break;
       end
       dat(:,:,j) = TifLink.read();
    end
    TifLink.close();
end