%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of Environmental Contours of Extreme Sea States
% Developed by Aubrey Eckert-Gallup, Cedric Sallaberry, and Ann Dallman
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1 - July 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Required Functions and Scripts:
%   princomp_rotation (Function)
%   dataorg (Function)
%   Comp2_bins (Function)
%   sigma_fits (Function)
%   mu_fits (Function)
%   iform (Function)
%   princomp_inv (Function)
%   steepness (Function)
%   dispersion_solver_NR_method (Function)
%   density_calc (Function)
%   Extreme_SeaState_Plots (Script)
%
% Variable Descriptions:
%   code_version= Current version of the code, used in the generation of
%                 output file names.
%   run_date    = Current date, used in the generation of output file
%                 names.
%   dataset_name= Name of data set to be considered, must be in .mat
%                 format and must include three vectors of input values: 
%                 DateNum, Hs, and Te or Tp.
%   depth       = Water depth at buoy under analysis.
%   size_bin    = Bin size that will be used to split Component 2 into
%                 bins after ranking according to values of Component 1.
%   nb_steps    = Discretization of the circle in the normal space used
%                 for inverse FORM calculation.
%   Time_SS     = Sea state duration (hours) of measurements in input. 
%   Time_r      = Desired return period (years) for calculation of
%                 environmental contour, can be a scalar or a vector.
%   SteepMax    = Optional user input of wave breaking steepness used to
%                 calculate steepness curves and adjust the extreme sea
%                 state contour accordingly.
%   DateNum     = Vector of timestamps for each measurement in the input.
%   Hs          = Vector of Hs values for each measurement in the input.
%   T           = Vector of Te or Tp values for each measurement in the 
%                 input. 
%	Comp1_Comp2 = Matrix with first column containing values of 
%                 component 1 and second column containing values of 
%                 component 2 after principal component rotation is 
%	              applied to Hs and T data.
%	coeff       = Principal component coefficients.
%	shift       = Shift applied to Component 2 to ensure that there aren't
%	              any negative values.
%	n_data      = Counts the number of measurements.
%   Rank_Comp1_Comp2 = Matrix with the rank(based on Component 1 value),
%                 Component 1, and corresponding Component 2 and DateNum
%                 measurements.
%   edges1      = Vector with the index of the edges of each bin of 
%                 size_bin up to the last bin.
%   histnum1    = Vector with number of values in each bin.
%   in_bin1     = Vector assigning a bin number to each measurement.
%   distnm_Comp1= Distribution chosen to fit the entire CDF of Component 1.
%   Comp1_pd    = Probability distribution object containing the fitted 
%                 Component 1 CDF.  
%   step_size   = Increment for CDF quantile calculation.
%   Comp1_ecdf  = CDF of Component 1 data.
%   Comp1_freq  = Matrix in which each column corresponds to a single bin
%                 and contains the values of Component 1 in each bin. 
%                 Note: The last column of this matrix may contain zeros  
%                 for rows that do not have measurements if the number of 
%                 bins is not a divisor of the number of measurements. 
%   Comp2_freq  = Matrix in which each column corresponds to a single bin
%                 and contains the values of Component 2 in each bin. 
%                 See note above.
%   Comp1_mean  = Vector of mean values of Component 1 for each bin.
%   distnm_Comp2= Distribution chosen to fit the CDF of Component 2 for 
%                 each bin.
%   mu_vals     = Mean of Component 2 for each bin based on fitted 
%                 distribution.
%   sigma_vals  = Standard deviation of Component 2 for each bin based on
%                 fitted distribution.
%   pd_names    = Name of probability distribution object for each bin of
%                 Component 2 so that each of these can be easily saved in 
%                 a loop over all bins.
%   Comp2_bin_pds = Structure containing the probability distribution 
%                 objects created by fitting the CDFs for Component 2 for 
%                 each bin.
%   sigma_param = Parameters of the function used to fit sigma as a
%                 function of the mean value of Component 1 for each bin.
%   sigma_fit   = Fitted value of sigma calculated at each mean value of
%                 Component 1 for each bin.
%   sigma_fcn   = Fitting function for sigma as a function of the mean
%                 value of Component 1 for each bin.
%   mu_param    = Parameters of the function used to fit mu as a
%                 function of the mean value of Component 1 for each bin.
%   mu_fit      = Fitted value of mu calculated at each mean value of
%                 Component 1 for each bin.
%   mu_fcn      = Fitting function for mu as a function of the mean
%                 value of Component 1 for each bin.
%   Comp1_R     = Calculated Component 1 values along the contour
%                 boundary.
%   Comp2_R     = Calculated Component 2 values along the contour
%                 boundary.
%   mu_R        = Calculated mean for the Component 2 values along the 
%                 contour boundary.
%   sigma_R     = Calculated standard deviation for the Component 2 values 
%                 along the contour boundary.
%   Hs_R        = Calculated Hs values along the contour boundary 
%                 following return to original input orientation.
%   T_R         = Calculated T values along the contour boundary 
%                 following return to original input orientation.
%   SteepH      = Wave height values calculated along the wave breaking
%                 curve.
%   T_vals      = Values of T at which the steepness values are
%                 calculated.
%   Hs_R_2      = Values of Hs along the contour boundary and along 
%                 the wave breaking steepness curve at values where the Hs 
%                 value is greater than the breaking steepness.
%   density_data_name = Name of file in which data density calculation will
%                 be stored.
%   Count_val   = Number of data points inside of a fixed neighborhood 
%                 around each point, used to describe the density of the 
%                 data for each point.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation setup and User Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

code_version = {'v_1'}; % Enter code version
run_date = {'07_21_15'};  % Enter date

dataset_name = {'NDBC46022_1996_2012_Hs_Te.mat'};
depth = 675; % Enter depth at measurement point (m)

size_bin = 150;  % Enter chosen bin size
nb_steps = 1000; % Enter discretization of the circle in the normal space  
                 % used for inverse FORM calculation
Time_SS = 1;   % Enter seastate duration (hours)
Time_r = [100];  % Enter return period (years) as a scalar or a row 
                       % vector.
SteepMax = 0.07; % Optional: enter estimate of breaking steepness; comment 
                 % this line to skip this step

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing NaN data, assigning T label depending on input (Te or Tp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(dataset_name{1}) % Loads timeseries of Hs, T, DateNum

tic % Start calculation timer

% Depending on input data for T (Tp or Te), assign data and plot labels:
if exist('Te','var')
    T = Te;
    TPlotLabel = 'Energy Period   T_e   [s]';
else if exist('Tp','var')
    T = Tp;
    TPlotLabel = 'Peak Period   T_p   [s]';
    end
end
    
% Remove NaN's:
NaNrem = find(isnan(Hs)|isnan(T)); % Find NaN data in either Hs or T
DateNum(NaNrem) = []; % Remove any NaN data from DateNum
Hs(NaNrem) = []; % Remove any NaN data from Hs
T(NaNrem) = []; % Remove any NaN data from T

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal Components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Applying principal component rotation')

[Comp1_Comp2,coeff,shift] = princomp_rotation(Hs,T);
% This function applies a rotation to the data using the method of
% principal components.

%%                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data organization and processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

[n_data,Rank_Comp1_Comp2,edges1,histnum1,in_bin1] = ...
    dataorg (Comp1_Comp2,DateNum,size_bin);
% This function organizes the data for new_extreme_wave_events.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting Component 1 distribution for whole data set 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Starting fitting distribition for Component 1')

distnm_Comp1 = 'InverseGaussian'; 
% Distribution chosen to fit entire set of Component 1. Options for this 
% choice can be found in the Matlab documentation for the inputs arguments 
% of the 'fitdist' function.

Comp1_pd = fitdist(Rank_Comp1_Comp2(:,2),distnm_Comp1);
% This function creates a probability distribution object by fitting a
% distribution to Component 1.

step_size = 1/(length(Rank_Comp1_Comp2)+1);
Comp1_ecdf = [step_size:step_size:(length(Rank_Comp1_Comp2)*step_size)]'; 
% Calculates the CDF of Component 1.


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Splitting Component 2 into bins of size size_bins according to sorted 
% Component 1 values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Splitting Component 2 into bins of size '...
        num2str(size_bin) ' according to Component 1'])

[Comp1_freq,Comp2_freq,Comp1_mean] = ...
    Comp2_bins(size_bin,edges1,histnum1,in_bin1,Rank_Comp1_Comp2);
% This function creates Component 2 bins and calculates related values. 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting each bin of Component 2 with a normal distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distnm_Comp2 = 'Normal';
display('Fitting normal distribution to each H_s bin')
mu_vals = zeros((length(edges1)-1),1);
sigma_vals = zeros((length(edges1)-1),1);
pd_names = cell((length(edges1)-1),1);
% Create probability distribution object for each bin by fitting binned
% values of Component 2 (Hs) with chosen distribution. Options for this 
% choice can be found in the Matlab documentation for the inputs arguments 
% of the 'fitdist' function.
for j = 1:(length(edges1)-1)
    pd_names{j} = strcat('Comp2_bin_pd',num2str(j));
    Index = find(Comp2_freq(:,j));
    Comp2_bins_pds.(pd_names{j})=fitdist(Comp2_freq(Index,j),distnm_Comp2); 
    % Save all probability distribution objects in a data structure
    mu_vals(j) = mean(Comp2_bins_pds.(pd_names{j})); % Calculation of mean (mu)
    sigma_vals(j) = std(Comp2_bins_pds.(pd_names{j})); % Calculation of standard 
                                                  % deviation (sigma)
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting mu and sigma as functions of mean Component 1 value for
% each bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Fitting distributions for mean and standard deviation')

[sigma_param,sigma_fit,sigma_fcn]= sigma_fits(Comp1_mean,sigma_vals);
% This function returns the fitting parameters for sigma as a function of 
% the mean of Component 1 for each bin fit by a linear approximation. 

[mu_param,mu_fit,mu_fcn] = mu_fits(Comp1_mean,mu_vals);
% This function returns the fitting parameters for mu as a function of 
% the mean of Component 1 for each bin fit by a linear approximation.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform the inverse FORM calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Performing inverse FORM calculation')

[Comp1_R,Comp2_R,mu_R,sigma_R] = ...
    iform(Time_r,Time_SS,nb_steps,Comp1_pd,mu_fcn,sigma_fcn);
% This function applies the inverse FORM technique to the data in order to
% estimate the extreme sea state contour for a given return period.


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return to original orientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Returning to original orientation')

[Hs_R,T_R] = princomp_inv(Comp2_R,Comp1_R,coeff,shift);

Hs_R = max(0,Hs_R); % Remove non-physical negative values for Hs from the 
                    % contour and replace them with zero


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate steepness curve -- if input is provided by user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('SteepMax','var')
    
display ('Calculating wave breaking steepness')

% Calculate the steepness curve for the given max steepness:
[SteepH,T_vals] = steepness(depth,max(T_R),SteepMax,[]);

% Calculate a second contour that includes the wave breaking steepness 
% curve at values where the Hs value on the contour is greater than the 
% breaking steepness
[SteepH_R] = steepness(depth,max(T_R),SteepMax,T_R);
Hs_R_2 = zeros(length(T_R),length(Time_r));
for j = 1:size(Time_r,2)
    for i = 1:length(T_R)
        if SteepH_R(i)<Hs_R(i,j)
            Hs_R_2(i,j) = SteepH_R(i);
        else
            Hs_R_2(i,j) = Hs_R(i,j);
        end
    end
end

else
    SteepH = [];
    Hs_R_2 = [];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Density Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

density_data_name = ...
    strcat(dataset_name{1}(1:(end-4)),'_countval_density.mat');

% Calculation of data density if this calculation hasn't already been made.
if isequal(exist(density_data_name,'file'),2)
    load(density_data_name);
else
    display('Density data does not exist in the current directory')
    display('Running density calculation')
 
    Count_val = density_calc(Hs,T,dataset_name);
end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting and Saving Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Plotting results')

run 'Extreme_SeaState_Plots.m'

display ('Saving results')

ExtremeSeaState_results = struct('Hs',Hs,'T',T,'DateNum',DateNum,...
    'Rank_Comp1_Comp2',Rank_Comp1_Comp2, 'Comp1_ecdf',Comp1_ecdf,...
    'Comp1_pd',Comp1_pd,'Comp1_freq',Comp1_freq, 'Comp2_freq',Comp2_freq,...
    'Comp2_bins_pds',Comp2_bins_pds,'histnum1',histnum1,'edges1',edges1,... 
    'Comp1_mean',Comp1_mean,'mu',mu_vals,'sigma',sigma_vals,'mu_fit',mu_fit,...
    'sigma_fit',sigma_fit,'mu_param',mu_param,'sigma_param',sigma_param,...
    'Comp1_R',Comp1_R,'Comp2_R',Comp2_R,'Hs_R',Hs_R,'Hs_R_2',Hs_R_2,...
    'T_R',T_R,'SteepH',SteepH,'depth',depth,...
    'princip_coeff',coeff,'shift',shift,'Time_r',Time_r);

save_name = strcat(dataset_name{1}(1:(end-4)),'_',code_version,'_',...
    run_date{1});
save(strcat(save_name{1},'.mat'),'ExtremeSeaState_results');

mkdir(save_name{1});
current_folder = pwd;

for i = 1:13;
    baseFileName = sprintf('figure_%d',i);
    figure(i);
    fullFileName = fullfile(strcat(current_folder,'\',save_name{1}),...
        baseFileName);
    saveas(gcf,fullFileName,'fig');
end

display('Calculation is complete')
toc;
elapsedTime_min = toc/60;
disp(['Elapsed time is ' num2str(elapsedTime_min) ' minutes'])


% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.