function pixels = linePixels(x1, y1, x2, y2)
    % Initialize the output array
    pixels = [];

    % Calculate differences
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);

    % Determine the number of steps needed
    steps = max(dx, dy)*8;
    
    % Calculate increment in x & y for each step
    Xinc = (x2 - x1) / steps;
    Yinc = (y2 - y1) / steps;
    
    % Initialize starting point
    X = x1;
    Y = y1;

    for i = 0:steps
        % Add the current pixel (using floor to determine the pixel it falls into)
        pixels = [pixels; round(X), round(Y)];
        
        % Move to the next point
        X = X + Xinc;
        Y = Y + Yinc;
    end
    
    % Remove duplicate pixels
    pixels = unique(pixels, 'rows');
end