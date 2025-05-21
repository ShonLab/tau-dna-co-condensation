function [y,n] = closest(A,x) % find x in A
[y,n] = deal(zeros(numel(x),1));
for i = 1:numel(x)
    n(i) = find(abs(A-x(i)) == min(abs(A-x(i))),1,'first');
    y(i) = A(n(i)); 
end