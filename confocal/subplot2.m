function h = subplot2(m,n,p,varargin)
    c = mod(p-1,n);
    r = floor((p-1)/n)+1;

    tmp = get(gcf,'DefaultAxesUnits');
    set(gcf,'DefaultAxesUnits','norm');
    h = axes('position',[c*(.99/n),1-r*(1/m),.99/n,.95/m],varargin{:}); 
    set(gcf,'DefaultAxesUnits',tmp);
end