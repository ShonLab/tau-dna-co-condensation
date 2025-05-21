function output = msum(input, dim)

if nargin == 1
    output = sum(input(:));
else
    output = input;
    for n = 1:numel(dim)
        output = sum(output, dim(n)); 
    end
end