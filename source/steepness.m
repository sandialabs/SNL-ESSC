function [SteepH,T_valsF] = steepness(depth,T_max,SteepMax,T_vals)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates a steepness curve to be plotted on an H vs T
% diagram.  First, the function calculates the wavelength based on the
% depth and T. The T vector can be the input data vector, or will be
% created below to cover the span of possible T values.
%
%   Syntax: [Steepness,T_mesh] = my_steepness(depth,T_max,Hs_max,T_vals,SteepMax)
%            where T_vals is optional input
%   Where: depth    = water depth [m]
%          T_max    = max value of period [sec]
%          SteepMax = wave breaking steepness (e.g., 0.07)
%          T_vals   = vector of T values [sec] (optional)
% 
%          T_valsF  = T values over which the steepness curve is defined [sec] 
%          SteepH   = the H values [m] that correspond to the T_mesh values 
%                   creating the steepness curve 
%
%
% Author Ann Dallman                                         Date 9/15/2014
%       based on many previous codes available online (mathworks exchange)
%       and one by Diana Bull
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Create steepness curves:

% If a vector is not given, create a finely spaced vector over which to
% calculate the H value at a given steepness:
if isempty(T_vals)
    T_vals = [0:.1:T_max]';
end
T_valsF = T_vals;

% Calculate the wavelength at a given depth at each value of T
lambdaT = NaN(size(T_vals));
for i = 1:length(T_vals(:,1))
    [~,lambdaT(i),~,~] = dispersion_solver_NR_method(depth,T_vals(i));
end

% Calculate the H values over which the steepness curve is defined:
SteepH = SteepMax.*lambdaT;
end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.