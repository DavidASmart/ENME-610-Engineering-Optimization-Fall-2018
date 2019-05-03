
function [X0] = genX0(d_L, d_U, t_L, t_U, L_L, L_U, N, Nd, Nt, NL, m)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Generates a set of start points for multi-start
%
%% Inputs:
%   d_L        lower bound for AUV inner diameter            (m)
%   d_U        upper bound for AUV inner diameter            (m)
%   t_L        lower bound for AUV hull thickness            (m)
%   t_U        upper bound for AUV hull thickness            (m)
%   L_L        lower bound for AUV length                    (m)
%   L_U        upper bound for AUV length                    (m)
%   N          number of points to generate
%   Nd         number of points in d (for grid only)
%   Nt         number of points in t (for grid only)
%   NL         number of points in L (for grid only)
%   m          method of point generation
%               1 - gaussian random
%               2 - halton psudo-random
%               3 - grid
%               4 - sphere
%
%% Outputs:
%   X0          set of start points
%

%% generate starting points
rng(1)
if m == 1
    X0 = makePaternRND(d_L, d_U, t_L, t_U, L_L, L_U, N);
    
elseif m == 2
    X0 = makePaternRND_Halton(d_L, d_U, t_L, t_U, L_L, L_U, N);
    
elseif m == 3
    X0 = makePaternGrid(d_L, d_U, t_L, t_U, L_L, L_U, Nd, Nt, NL);
    
else
    X0 =  makePaternS(d_L, d_U, t_L, t_U, L_L, L_U, N);
end

end
