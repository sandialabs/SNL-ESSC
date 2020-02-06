function [Comp1_Comp2,coeff,shift] = princomp_rotation(Hs,T)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function applies a rotation to the data using the method of
% principal components.
%
%   Syntax: [Comp1_Comp2,coeff,shift] = princomp_rotation(Hs,T)
%   Variables:
%   Hs          = Vector of Hs values for each measurement in the input.
%   T           = Vector of T values for each measurement in the input.
%   Hs_T_cat    = Concatenation of Hs, T vectors into one matrix with two
%                 columns, required for the use of Matlab function
%                 princomp.
%   coeff       = Principal component coefficients.
%   Comp1_Comp2 = Matrix with first column containing values of component 1
%                 and second column containing values of component 2 after 
%                 principal component rotation is applied to Hs and T
%                 data.
%   shift       = Shift applied to Component 2 to ensure that there aren't
%                 any negative values.
%
% Author: Aubrey Eckert-Gallup
% Date: 01/13/14
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hs_T_cat = horzcat (Hs,T);

% Calculate coefficients needed to rotate data using method of principal
% components
coeff = princomp(Hs_T_cat);

% Apply principal component roation to data
Comp1_Comp2 = Hs_T_cat * coeff;

shift = 0; % Preallocate shift

% Apply a shift to the data to ensure that there aren't any negative values
% of Component 2 (important for later calculations)
if any(Comp1_Comp2(:,2)<0)
    shift = abs(min(Comp1_Comp2(:,2)))+0.1;
    Comp1_Comp2(:,2) = shift+Comp1_Comp2(:,2);
end
end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.