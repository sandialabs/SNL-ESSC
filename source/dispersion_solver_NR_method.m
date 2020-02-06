function [k,lambda,Cp,Cg] = dispersion_solver_NR_method(h,T)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves the dispersion relation for water waves using the
% Newton-Raphson method. All outputs are solved for exactly using: 
% (w^2*h/g=kh*tanh(kh))
%
% Approximations that could be used in place of this code for deep and 
% shallow water, as appropriate:
% deep water: h/lambda >= 1/2, tanh(kh)~1, lambda = (g.*T.^2)./(2*.pi)
% shallow water: h/lambda <= 1/20, tanh(kh)~kh, lambda = T.*(g.*h)^0.5
%
%   Syntax: [k,lambda,C,Cg] = dispersion_solver_NR_method(h,T)
%   Where: h = water depth [m]
%          T = period [sec]
%          k = wave number [1/m]
%          lambda = wavelength [m]
%          Cp = phase velocity [m/s]
%             = w/k = (g/k*tanh(kh))^.5 = lambda/T
%          Cg = group velocity [m/s] 
%             = del(omega)/del(k) = nCp; n = 0.5(1+2kh/sinh(2kh))
%          kh = k*h 
%          w = omega = 2*pi / T = angular frequency [rad/s] 
%
%
% Author Ann Dallman                                         Date 9/15/2014
%       based on many previous codes available online (mathworks exchange)
%       and one by Kelley Ruehl
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


g = 9.81;  % [m/s^2]
omega = ((2*pi)./T); % [rad/sec]

%Pre-allocate vectors:
k = NaN(length(T),1);
lambda = NaN(length(T),1);
Cp = NaN(length(T),1);
n = NaN(length(T),1);
Cg = NaN(length(T),1);

% Calculate over the timeseries (vector of T)
for j = 1:length(T)
    
    %Initialize kh using Eckert 1952 (mentioned in Holthuijsen pg. 124)
    kh = omega(j).^2 .* h ./ (g * (tanh(omega(j).^2.*h./g)).^0.5 );
    
    %%Find solution using the Newton-Raphson Method.
    for i = 1:1000
        kh0 = kh; 
        f0 = omega(j).^2.*h./g - kh0.*tanh(kh0); %function w^2*h/g=kh*tanh(kh) --> f0 = w^2*h/g - kh*tanh(kh)
        df0 = -tanh(kh0) - kh0*(1-tanh(kh0)^2); %first derivative with respect to kh
        kh = -f0/df0 + kh0;
        f = omega(j).^2.*h./g - kh.*tanh(kh);
        if abs(f0-f) < 10^(-6), break, end
    end
    
    %Output
    k(j) = kh/h;
    lambda(j) = (2*pi)/k(j);
    Cp(j) = lambda(j)./T(j); %also, Cp = w/k = (g/k*tanh(kh))^.5 Holthuijsen pg. 125
    n(j) = 0.5*(1+((2*k(j)*h)/sinh(2*k(j)*h))); % Holthuijsen pg. 127
    Cg(j) = n(j).*Cp(j);
    
    clear kh i kh0 f0 df0 f
end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.
