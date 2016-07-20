function out = goodmapping(in,FOV)
% out = goodmapping(in, FOV)
%
% Compute 1D to 2D using the the good thomas mapping. Assumes the 2D dimensions are
% co-prime
%
% Inputs: 
%     in      -   Input 1D data
%     FOV     -   Field of view (2D dimension)
%
% Outputs:
%     out     -  Output 2D data
%       
% (c) Frank Ong 2013

[x,y] = ndgrid(1:FOV(2):length(in),1:FOV(1):length(in));

ind = mod(x+y-2,length(in))+1;

out = in(ind);