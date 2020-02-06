function [original1,original2] = princomp_inv(princip_data1,princip_data2,coeff,shift)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the inverse of the principal component rotation in
% order to return to the original orientation of the data.
%
%   Syntax: [original1,original2] = princomp_inv(princip_data,coeff,shift)
%	Variables:
%	princip_data = Matrix with first column containing values of 
%                  component 1 and second column containing values of 
%                  component 2 after principal component rotation is 
%	               applied to Hs and T data.
%	coeff        = Principal component coefficients.
%	shift        = Shift applied to Component 2 to ensure that there aren't
%	               any negative values.
%   original1    = Values of princip_data column 1 rotated back into 
%                  original input space.
%   original2    = Values of princip_data column 2 rotated back into 
%                  original input space.
% Author: Aubrey Eckert-Gallup
% Date: 01/13/14
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
original1 = zeros(size(princip_data1,1),size(princip_data1,2));
original2 = zeros(size(princip_data1,1),size(princip_data1,2));

for i = 1:size(princip_data1,2)
original1(:,i) = ((coeff(1,2).*(princip_data1(:,i)-shift))+(coeff(1,1).*princip_data2(:,i)))./(coeff(1,2)^2+coeff(1,1)^2);
original2(:,i) = ((coeff(1,2).*princip_data2(:,i))-(coeff(1,1).*(princip_data1(:,i)-shift)))./(coeff(1,2)^2+coeff(1,1)^2);
end

end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.