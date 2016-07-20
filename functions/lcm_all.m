function out = lcm_all(in)
% out = lcm_all(in)
%
% Compute least common multiplier for all elements in input. Matlab lcm only
% does pair-wise
%
% Inputs: 
%     in      -   Input
%
% Outputs:
%     out     -  lcm(in)
%       
% (c) Frank Ong 2013


out = in(1);
for i = 2:length(in)
    out = lcm(out,in(i));
end

