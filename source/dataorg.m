function [n_data,Rank_Comp1_Comp2,edges1,histnum1,in_bin1] = dataorg (Comp1_Comp2,DateNum,size_bin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function organizes the data for new_extreme_wave_events and 
% determines the bins of size "size_bin" for later use.
%
%   Syntax: [n_data,Rank_Comp1_Comp2,edges1,histnum1,in_bin1] = 
%            dataorg (Comp1_Comp2,DateNum,size_bin)
%	Variables:
%   Comp1_Comp2 = Matrix with first column containing values of 
%                 component 1 and second column containing values of 
%                 component 2 after principal component rotation is 
%	              applied to Hs and Te data.
%   DateNum     = Vector of timestamps for each measurement in the input.
%   size_bin    = Bin size that will be used to split Component 2 into
%                 bins after ranking according to values of Component 1.
%	n_data      = Counts the number of measurements.
%   Comp1_Comp2_DateNum = Concatenated matrix with columns of Component 1,
%                 Component 2, and DateNum.
%   sorted_by_Comp1 = Sorts the Comp1_Comp2_DateNum matrix by the first
%                 column (values of Component 1) in ascending order.
%   rank_val    = Rank of each measurement from 1 to the number of
%                 measurements.
%   Rank_Comp1_Comp2 = Matrix with the rank(based on component 1), values 
%                 of Component 1 sorted in ascending order, and 
%                 corresponding Component 2 and DateNum measurements.
%   edges1      = Vector with the index of the edges of each bin of 
%                 size_bin up to the last bin.
%   histnum1    = Vector with number of values in each bin.
%   in_bin1     = Vector assigning a bin number to each measurement.
% 
% Author: Aubrey Eckert-Gallup
% Date: 01/13/14
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_data = length(Comp1_Comp2); % Count data
Comp1_Comp2_DateNum = cat(2,Comp1_Comp2,DateNum); % Create Component1, Component2, DateNum matrix

rank_val = [1:1:n_data]'; % Create matrix of integers up to number of 
                          % data points corresponding to the rank
Rank_Comp1_Comp2 = cat(2,rank_val,sortrows(Comp1_Comp2_DateNum,1)); % Creates matrix including rank 
                                           % in first column and values
                                           % sorted by Component 1
edges1 = [1:size_bin:size_bin*ceil(n_data/size_bin),n_data+1]; % Creates edges of bins of size of 
                                     % size_bin up to last bin

[histnum1,in_bin1] = histc(rank_val,edges1); % Creates a matrix with number 
% of values in each bin and a matrix assigning bin number to each value index
end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.