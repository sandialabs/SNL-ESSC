function [mu_param,mu_fit,mu_fcn] = mu_fits(Comp1_mean,mu_vals)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the fitting parameters for mu as a function of 
% the mean of Component 1 for each bin fit by a linear approximation.
%
%   Syntax: [mu_param,mu_fit,mu_fcn] = mu_fits(Comp1_mean,mu)
%	Variables:
%   Comp1_mean  = Vector of mean values of Component 1 for each bin.
%   mu          = Mean of Component 2 for each bin based on fitted 
%                 distribution.
%   mu_param    = Parameters of the function used to fit mu as a
%                 function of the mean value of Component 1 for each bin.
%   mu_fit      = Fitted value of mu calculated at each mean value of
%                 Component 1 for each bin.
%   mu_fcn      = Fitting function for mu as a function of the mean
%                 value of Component 1 for each bin.
%   xData       = Comp1_mean values checked for errors before fitting 
%                 proceeds.
%   yData       = Mu values checked for errors before fitting proceeds.
%   ft          = Fit type.
%   opts        = Fitting options.
%   fitresult   = Fit object created by fitting the given function to the
%                 data.
%
% Author: Aubrey Eckert-Gallup
% Date: 01/13/14
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xData, yData] = prepareCurveData( Comp1_mean, mu_vals );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf -Inf];

% Fit model to data.
fitresult = fit( xData, yData, ft,opts);
mu_param = coeffvalues(fitresult);

mu_fcn = @(x) mu_param(1).*x + mu_param(2);
mu_fit = mu_fcn(Comp1_mean);
end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.