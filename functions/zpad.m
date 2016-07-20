function res = zpad(x,s)
%
%  res = zpad(x,[sx,sy,sz,st])
%  same as the previous example
%
%
% (c) Michael Lustig 2007

if nargin < 2
    error('must have a target size')
end

m = size(x);
if length(s) < length(m)
    s = [s, ones(1,length(m)-length(s))];
end
if length(s) > length(m)
    m = [m, ones(1,length(s)-length(m))];
end

if sum(m==s)==length(m)
    res = x;
    return;
end

res = zeros(s);

for n=1:length(s)
    idx{n} = floor(s(n)/2)+1+ceil(-m(n)/2) : floor(s(n)/2)+ceil(m(n)/2);
end
% this is a dirty ugly trick
cmd = 'res(idx{1}';
for n=2:length(s)
    cmd = sprintf('%s,idx{%d}',cmd,n);
end
cmd = sprintf('%s)=x;',cmd);
eval(cmd);





