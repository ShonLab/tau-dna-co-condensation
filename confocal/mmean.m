function output = mmean(input,dim)

if nargin == 1
    output = mean(input(:));
else
    for n = 1:numel(dim)
        input = mean(input, dim(n)); 
    end
    output = input;
end