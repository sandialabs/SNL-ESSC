function [Count_val] = density_calc(Hs,Te,dataset_name)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the number of points in a neighborhood around
% each point to create a representation of data density.
%
%   Syntax: [Count_val] = density_calc(Hs,Te,dataset_name)
%   Variables:
%   Hs          = Vector of Hs values for each measurement in the input.
%   Te          = Vector of Te values for each measurement in the input. 
%   dataset_name= Name of data set to be considered, must be in .mat
%                 format and must include three vectors of input values: 
%                 DateNum, Hs, Te.
%   mean_Hs     = Mean of Hs.
%   mean_Te     = Mean of Te.
%   stdv_Hs     = Standard deviation of Hs.
%   stdv_Te     = Standard deviation of Te.
%   Hs_2        = Values of Hs normalized by subtracting the mean and
%                 dividing by the standard deviation.
%   Te_2        = Values of Te normalized by subtracting the mean and
%                 dividing by the standard deviation.
%   Count_val   = Number of data points inside of a fixed neighborhood 
%                 around each point, used to describe the density of the 
%                 data for each point.
%   dist        = Radius of the neighborhood around each point.
%   disp_val    = Value used to display percent completion due to long
%                 calculation time.
%   save_name   = Name of file used to save calculated density data
%
% Author: Aubrey Eckert-Gallup
% Date: 01/13/14
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_Hs = mean(Hs); % Calculation of mean Hs value
mean_Te = mean(Te); % Calculation of mean Te value
stdv_Hs = std(Hs); % Calculation of Hs standard deviation
stdv_Te = std(Te); % Calculation of Te standard deviation

Hs_2 = zeros(length(Hs),1); 
Te_2 = zeros(length(Hs),1);

% Normalization
for i = 1:length(Hs);
    Hs_2(i) = (Hs(i)-mean_Hs)/stdv_Hs; 
    Te_2(i) = (Te(i)-mean_Te)/stdv_Te;    
end 

Count_val = zeros(length(Hs),1);

% Count values inside of elipsoid for each point
dist = 2000/length(Hs);
display('Starting count of neighborhood values')
disp_val = floor(length(Hs)/10);
for i = 1:length(Hs)
    for j = 1:length(Hs)
        if sqrt((Hs_2(i)-Hs_2(j))^2+(Te_2(i)-Te_2(j))^2)<dist
            Count_val(i) = Count_val(i)+1;
        end
    end
    if i == disp_val
        display(strcat(num2str(ceil(disp_val/length(Hs)*100)),'% complete'))
        disp_val = disp_val + floor(length(Hs)/10);
    end
end

save_name = strcat(dataset_name{1}(1:(end-4)),'_countval_density.mat');
save(save_name,'Count_val');
end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.