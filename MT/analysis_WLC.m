%% load data
% set parameters
clear; clc;
Rbead = 1400; % MB radius in nm
T = 300; % temperature
pixel_size = 80; % in nm
correction_factor = 0.878; % correction factor for refractive index mismatch (n_water/n_oil = 0.878)
% F0 = 0.05759; z0 = 20.94; A1 = 91.75; d1 = 0.9161; A2 =  53.88; d2 = 2.647; %JOVE
%  F0 = 0.1; z0 = 20.6; A1 = 194.2; d1 = 1.009; A2 = 68.45; d2 =2.779; % M270 (20221214 for POSTECH fresh setup,TY calibrated with N52 magnets;)
% F0 = 0.01819; d1 = 3.3; d2 = 1.573; A1 = 22.23; A2 = 50; z0 = 19.83; 
% F0 = 0.03448; d1 = 1.106; d2 = 2.615; A1 = 35.56; A2 = 43.14; z0 = 50; % 2023 MIE calib. JJ
F0 = 1e-6; d1 = 1.361; d2 = 3.117; A1 = 54.46; A2 = 27.35; z0 = 49.55; % 20240120. JJ

% read data
folders = dir('WLC\raw*');
pths = arrayfun(@(i) {['WLC\',folders(i).name,'/']}, 1:numel(folders)); npth = numel(pths); %pths =raw*n
[finfo,fname,nbead,nframe,fps,Roff,ori,f,t,t2,d,R,P,rz,rx,ry,x,y,z,dx,dy,dz] = deal(cell(npth,1));
nfile = zeros(npth,1);
for p = 1:npth
    
    disp(p);
    finfo{p} = dir([pths{p},'r*.xls']); % data files
    nfile(p) = numel(finfo{p});
    [fname{p},Roff{p},ori{p},f{p},d{p},R{p},P{p},x{p},y{p},z{p},dx{p},dy{p},dz{p}] = deal(cell(nfile(p),1));
    [nbead{p},fps{p},nframe{p}] = deal(zeros(nfile(p),1));
    for n = 1:nfile(p)
        disp([int2str(n/nfile(p)*100),'% of ',pths{p}(1:end-1),'...']);
        fname{p}{n} = finfo{p}(n).name;
        
        % read calibration info (c###.xls)
        calinfo = dlmread([pths{p},'c',fname{p}{n}(2:4),'.fps']);
        fps{p}(n) = calinfo(1); % sampling rate in Hz (1200 here)
        Roff{p}{n} = calinfo(2,:); % tether offset (see: 10.1126/sciadv.aav1697)
        ori{p}{n} = calinfo(3,:); % +1/-1 for left-/right-tilted MB, respectively
        
        % read motor data (s###-###.xls)
        dat = dlmread([pths{p},regexprep(fname{p}{n},'r','s')]);
        t2{p}{n} = dat(:,1); % timestamp for motor data in seconds
        M{p}{n} = dat(:,2); % magnet distance in mm
        F{p}{n} = F0 + A1*exp(-(z0-M{p}{n})/d1) + A2*exp(-(z0-M{p}{n})/d2);
        R{p}{n} = (dat(:,3)-floor(dat(:,3)))*360; % magnet orientation in degrees
        P{p}{n} = dat(:,4); % piezo position in μm
        
        % read bead data (r###-###.xls)
        dat = dlmread([pths{p},fname{p}{n}]);
        nframe{p}(n) = size(dat,1); % # of data points
        f{p}{n} = dat(:,1); % frame number
        t{p}{n} = f{p}{n}/fps{p}(n); % timestamp for bead data in seconds
        dat = dat(:,2:end); % remove the first column
        dat(:,[1:3:end,2:3:end]) = dat(:,[1:3:end,2:3:end]) - repmat(mean(dat(31:60,[1:3:end,2:3:end]),1),[nframe{p}(n),1]); % subtract xy offset
        nbead{p}(n) = size(dat,2)/3-1; % # of MB (one is for RB)
        
        % bead position
        rx{p}{n} = dat(:,1)*pixel_size;
        ry{p}{n} = dat(:,2)*pixel_size;
        rz{p}{n} = dat(:,3);
        x{p}{n} = dat(:,4:3:end)*pixel_size;
        y{p}{n} = dat(:,5:3:end)*pixel_size;
        z{p}{n} = dat(:,6:3:end);
        
        % MB position relative to RB
        dx{p}{n} = (x{p}{n}-repmat(rx{p}{n},[1,nbead{p}(n)]));
        dy{p}{n} = (y{p}{n}-repmat(ry{p}{n},[1,nbead{p}(n)]));
        dz{p}{n} = (z{p}{n}-repmat(rz{p}{n},[1,nbead{p}(n)]))*correction_factor;
      
        % synchronize motor data to bead data
        F{p}{n} = interp1(t2{p}{n},F{p}{n},t{p}{n});
        R{p}{n} = interp1(t2{p}{n},R{p}{n},t{p}{n});
        P{p}{n} = interp1(t2{p}{n},P{p}{n},t{p}{n});       
    end
end
clear dat;
save('analysis_WLC');