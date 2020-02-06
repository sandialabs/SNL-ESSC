function [sigma_param,sigma_fit,sigma_fcn,rho]= sigma_fits(Comp1_mean,sigma_vals);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the fitting parameters for sigma as a function of 
% the mean of Component 1 for each bin fit by a linear approximation. 
%
%   Syntax: [sigma_param,sigma_fit,sigma_fcn]= sigma_fits(Comp1_mean,sigma)
%   Variables:
%   Comp1_mean  = Vector of mean values of Component 1 for each bin.
%   sigma       = Standard deviation of Component 2 for each bin based on
%                 fitted distribution.
%   sigma_param = Parameters of the function used to fit sigma as a
%                 function of the mean value of Component 1 for each bin.
%   sigma_fit   = Fitted value of sigma calculated at each mean value of
%                 Component 1 for each bin.
%   sigma_fcn   = Fitting function for sigma as a function of the mean
%                 value of Component 1 for each bin.
%   xData       = Comp1_mean values checked for errors before fitting
%                 proceeds.
%   yData       = Sigma values checked for errors before fitting proceeds
%   ft          = Fit type.
%   opts        = Fitting options.
%   fitresult   = Fit object created by fitting the given function to the
%                 data.
%
% Author: Aubrey Eckert-Gallup
% Date: 01/13/14
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xvals = Comp1_mean; 
yvals = sigma_vals;
s_0 = [0.1 0.1 0.1]; % Set intial guess

rho = 1; % Set intial penalty value
epsilon = 10^-5; % Set tolerance, very small values 
                 % (i.e., smaller than 10^-5) may cause instabilities

[Beta1,Beta2] = betafcn(s_0); % Set initial beta values using beta function
fun = @(x,s) s(1).*x.^2 + s(2).*x + s(3); % 2nd order polynomial function
objfun = @(s) sum((fun(xvals,s)-yvals).^2); % Sum of least squares 
                                            % objective function
objfun_penalty = @(S) objfun(S) + Beta1*(-S(3))^2 +Beta2*(-S(3)+(S(2)^2)/(4*S(1)))^2; 
% Objective function with penalty functions for both constraints

s_1 = fminsearch(objfun_penalty,s_0); % Initial search for minimum value using initial guess

% While either the difference between iterations or the difference in the
% objective function evaluation is greater than the tolerance, continue
% iterating
while any(abs(s_1-s_0)>epsilon) && abs(objfun(s_1)-objfun(s_0))>epsilon 
    s_0 = s_1;
    [Beta1,Beta2] = betafcn(s_0); % Set initial beta values
    objfun_penalty = @(S) objfun(S) + Beta1*(-S(3))^2 +Beta2*(-S(3)+(S(2)^2)/(4*S(1)))^2;
    s_1 = fminsearch(objfun_penalty,s_0);
    rho = 10*rho;
end
sigma_param = s_1;

function [Beta1,Beta2] = betafcn(s)
    if -s(3)<=0 % If first constraint is satisfied, do not apply a penalty
        Beta1 = 0;
    else % If not, apply penalty
        Beta1 = rho;
    end
    if -s(3)+(s(2)^2)/(4*s(1)) <= 0 % If second constraing is satisfied, do not apply a penalty
        Beta2 = 0;
    else % If not, apply penalty
        Beta2 = rho;
    end
end

sigma_fcn = @(x) sigma_param(1).*x.^2 + sigma_param(2).*x+sigma_param(3);
% Create sigma function object using the parameters that have been
% calculated
sigma_fit = sigma_fcn(Comp1_mean);
% Save the fitted values for plotting later

end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.