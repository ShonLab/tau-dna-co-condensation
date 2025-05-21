function hFig = maxfig(hFig,varargin)
pos_screen = get(0,'MonitorPositions');
[~,sel] = max(pos_screen(:,3));
pos_fig = pos_screen(sel,:);

if nargin == 0
    hFig = figure;
elseif nargin>=1
    for f = 1:numel(hFig)
        if ~ishandle(hFig(f))
            figure(hFig(f));
        end
    end
end

for f = 1:numel(hFig)
    set(hFig(f),'outerposition',pos_fig);
    drawnow % Required to avoid Java errors
    jFig = get(hFig(f), 'JavaFrame');
    jFig.setMaximized(true);
end
sfigure(hFig(end));

if nargin>=2
    set(hFig(:),varargin{:});
end
