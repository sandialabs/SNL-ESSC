function [Comp1_R,Comp2_R,mu_R,sigma_R] = iform(Time_r,Time_SS,nb_steps,Comp1_pd,mu_fcn,sigma_fcn)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function applies the inverse FORM technique to the data in order to
% estimate the extreme event contour for a given return period.
%
%   Syntax: [Comp1_R,Comp2_R,mu_R,sigma_R] = 
%            iform(Time_r,Time_SS,nb_steps,Comp1_pd,mu_fcn,sigma_fcn)
%	Variables:
%   Time_r      = Desired return period (years) for calculation of
%                 environmental contour.
%   Time_SS     = Sea state duration (hours) of measurements in input. 
%   nb_steps    = Discretization of the circle in the normal space used
%                 for inverse FORM calculation.
%   Comp1_pd    = Probability distribution object containing the fitted 
%                 Component 1 CDF.  
%   mu_fcn      = Fitting function for mu as a function of the mean
%                 value of Component 1 for each bin.
%   sigma_fcn   = Fitting function for sigma as a function of the mean
%                 value of Component 1 for each bin.
%   p_f         = Failure probability for the desired return period  
%                (Time_r) given the duration of the measurements (Time_SS).
%   beta        = Radius of reliability (distance from the most likely
%                 point in standard normal space).
%   theta       = Discretization of the circle based on the desired number
%                 of points.
%   U1          = Discrete values of X along the circumference of the 
%                 circle in standard normal space.  
%   U2          = Discrete values of Y along the circumference of the 
%                 circle in standard normal space.
%   Comp1_R     = Calculated Component 1 values along the extreme event
%                 boundary.
%   Comp2_R     = Calculated Component 2 values along the extreme event
%                 boundary.
%   mu_R        = Calculated mean for the Component 2 values along the 
%                 extreme event boundary.
%   sigma_R     = Calculated standard deviation for the Component 2 values 
%                 along the extreme event boundary.
%
% Author: Aubrey Eckert-Gallup
% Date: 01/13/14
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_f = 1 ./ (365.*(24./Time_SS).*Time_r);   % Failure probability
beta = norminv((1 - p_f),0,1);   % Reliability

% Vary U1,U2 along circle sqrt(U1^2+U2^2) = beta
theta = linspace(0,2*pi,nb_steps);

for i = 1:size(Time_r,2)
    U1(:,i) = beta(1,i)*cos(theta);
    U2(:,i) = beta(1,i)*sin(theta);
end

Comp1_R = icdf(Comp1_pd,normcdf(U1)); % Calculate the inverse CDF in order to
% estimate the Component 1 values along the extreme event boundary
 
% Use the approximation parameters previously calculated for
% mu and sigma to determine the values for mu and sigma of Component 2 along the
% extreme event contour.

% Apply mu fitting function
mu_R = mu_fcn(Comp1_R);
% Apply sigma fitting function
sigma_R = sigma_fcn(Comp1_R);
% Use the inverse of the normal distribution CDF to determine the Component 2 values
% for the extreme event contour.
Comp2_R = norminv(normcdf(U2),mu_R,sigma_R);

end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.