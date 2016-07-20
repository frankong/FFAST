function out = crtmapping(in)
% out = crtmapping(in)
%
% Compute 2D to 1D using the the crt mapping. Assumes the 2D dimensions are
% co-prime
%
% Inputs: 
%     in      -   Input 2D data
%
% Outputs:
%     out     -  Output 1D data
%       
% (c) Frank Ong 2013

out = zeros(1,numel(in));

for ind = 1:numel(in)
    i = mod ( ind-1 , size(in, 1) ) + 1;
    j = mod ( ind-1 , size(in, 2) ) + 1;
    out(ind) = in(i,j);
end