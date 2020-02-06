function [Comp1_freq,Comp2_freq,Comp1_mean] = Comp2_bins(size_bin,edges1,histnum1,in_bin1,Rank_Comp1_Comp2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates Component 2 bins and calculates related values.
%
%   Syntax: [Comp1_freq,Comp2_freq,Comp1_mean] = 
%            Comp2_bins(size_bin,edges1,histnum1,in_bin1,Rank_Comp1_Comp2)
%   Variables:
%   size_bin    = Bin size that will be used to split Component 2 into
%                  bins after ranking according to values of Component 1.
%   edges1      = Vector with the index of the edges of each bin of 
%                 size_bin up to the last bin.
%   histnum1    = Vector with number of values in each bin.
%   in_bin1     = Vector assigning a bin number to each measurement.
%   Rank_Comp1_Comp2 = Matrix with the rank(based on Component 1 value),
%                 Component 1, and corresponding Component 2 and DateNum
%                 measurements.
%   Index       = Index of values that are in a particular bin
%   Comp1_freq  = Matrix in which each column corresponds to a single bin
%                 and contains the values of Component 1 in each bin. 
%                 Note: The last column of this matrix may contain zeros 
%                 for rows that do not have measurements if the number of 
%                 bins is not a divisor of the number of measurements 
%   Comp2_freq  = Matrix in which each column corresponds to a single bin
%                 and contains the values of Component 2 in each bin. 
%                 See note above.
%   Comp1_mean  = Vector of mean values of Component 1 for each bin
%   
% Author: Aubrey Eckert-Gallup
% Date: 01/13/14
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Comp1_freq = zeros(size_bin,(length(edges1))-1); % Creates an empty matrix  
                                              % for storing bin values
Comp2_freq = zeros(size_bin,(length(edges1))-1); % Creates an empty matrix  
                                              % for storing bin values
Comp1_mean = zeros(length(edges1)-1,1); % Create matrix to store mean of 
                                        % Component 1 for all bins
                                     
% Calculate Component 2 CDF for each bin
for i = 1:(length(edges1)-1) %Loops over all bins
    Index = find(in_bin1 == i); % Finds the index of values that are in 
                                 % bin hk
    
    Comp1_freq(1:histnum1(i),i) = sortrows(Rank_Comp1_Comp2(Index,2));
    Comp2_freq(1:histnum1(i),i) = sortrows(Rank_Comp1_Comp2(Index,3));
    Comp1_mean(i) = sum(Comp1_freq(1:histnum1(i),i))/histnum1(i);
end
end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.