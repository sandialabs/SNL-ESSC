%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extreme Sea State Plots
% Aubrey Eckert-Gallup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script plots 13 figures for the Extreme Sea States code. These
% plots are ordered and organized into sections based on the timing of the 
% calculations made to create them. 

% Figure 1: Plotting the Component 1 distribution fit with original CDF.
figure(1)
set(gcf,'DefaultAxesFontSize',12);
box on
plot(Rank_Comp1_Comp2(:,2),Comp1_ecdf,Rank_Comp1_Comp2(:,2),cdf(Comp1_pd,Rank_Comp1_Comp2(:,2)),'r')
xlabel('Component 1');
ylabel('CDF: Prob(Component 1 \leq C_1)');
legend('Original CDF','Inverse Gaussian Fit','location','SouthEast');
grid on

% Figure 2: Plotting the Component 1 distribution fit with original CDF 
% zoomed in to the top of the distribution.
figure(2)
set(gcf,'DefaultAxesFontSize',12);
grid on
box on
hold on
plot(Rank_Comp1_Comp2(:,2),Comp1_ecdf,Rank_Comp1_Comp2(:,2),cdf(Comp1_pd,Rank_Comp1_Comp2(:,2)),'r')
xlabel('Component 1');
ylabel('CDF: Prob(Component 1 \leq C_1)');
legend('Original CDF','Inverse Gaussian Fit','Location','SouthEast');
axis([0 10 0.98 1])
axis 'auto x'
hold off

% Figure 3: Plotting the Component 1 distribution fit with original CDF 
% zoomed in to the bottom of the distribution.
figure(3)
set(gcf,'DefaultAxesFontSize',12);
grid on
box on
hold on
plot(Rank_Comp1_Comp2(:,2),Comp1_ecdf,Rank_Comp1_Comp2(:,2),cdf(Comp1_pd,Rank_Comp1_Comp2(:,2)),'r')
xlabel('Component 1');
ylabel('CDF: Prob(Component 1 \leq C_1)');
legend('Original CDF','Inverse Gaussian Fit','Location','SouthEast');
axis([0 2 0 0.2])
axis 'auto x'
hold off

%%
% Figure 4: Plotting the Component 2 CDF for each bin.
figure (4)
set(gcf,'DefaultAxesFontSize',12);
grid on
box on
hold on
for hk = 1:(length(edges1)-1)
    step_size = 1/((histnum1(hk))+1);
    quantiles = [step_size:step_size:((histnum1(hk))*step_size)];
    Index = find(Comp2_freq(:,hk));
    plot(Comp2_freq(Index,hk),quantiles)
    xlabel('Component 2');
    ylabel('CDF: Prob(Component 2 \leq C_2)');
end 
hold off

% Figure 5: Plotting the Component 2 CDF for each bin with distribution 
% fitting for five selected bins
figure (5)
set(gcf,'DefaultAxesFontSize',12);
grid on
box on
hold on
for hk = floor((length(edges1)-1)/5).*[1 2 3 4 5]
    step_size = 1/((histnum1(hk))+1);
    quantiles = [step_size:step_size:((histnum1(hk))*step_size)]';
    Index = find(Comp2_freq(:,hk));
    S = strcat('Comp2_bin_pd',num2str(hk));
    plot(Comp2_freq(Index,hk),quantiles,Comp2_freq(Index,hk),cdf(Comp2_bins_pds.(S),Comp2_freq(Index,hk)),'r');
    xlabel('Component 2');
    ylabel('CDF: Prob(Component 2 \leq C_2)');
    legend('Original CDF','Normal Fit','location','SouthEast');
end
hold off

%%
%Figure 6: Display mu and sigma fits
figure(6); set(gcf,'DefaultAxesFontSize',12);
subplot(2,1,1); plot(Comp1_mean,mu_vals,'k.',Comp1_mean,mu_fit,'k-');
xlabel('Component 1');
ylabel({'Mean, \mu'});
subplot(2,1,2); plot(Comp1_mean,sigma_vals,'k.',Comp1_mean,sigma_fit,'k-');
xlabel('Component 1');
ylabel({'Standard Deviation, \sigma'});

%%
% Figure 7: Plot the Component1 and Component 2 data and corresponding 
% extreme event boundary.
figure (7)
hold on;
box on
set(gcf,'DefaultAxesFontSize',12);
plot(Comp1_R,Comp2_R,'k');
plot(Comp1_Comp2(:,1), Comp1_Comp2(:,2), '.')
title('Extreme Sea State Contour within Principal Component Axis')
xlabel('Component 1');
ylabel('Component 2');
for i = 1:length(Time_r)
    [multmax,indmultmax] = max(Comp1_R(:,i).*Comp2_R(:,i));
    Comp2loc = Comp2_R(indmultmax,i);
    Comp1loc = Comp1_R(indmultmax,i);
    text(Comp1loc,Comp2loc,sprintf('%d yr',Time_r(i)),'FontWeight','bold','BackgroundColor','w','HorizontalAlignment','left','VerticalAlignment','middle','Rotation',-50);
end
%legend(sprintf('%d year contour',Time_r));
hold off

%%
% Figure 8: Plot the Hs and T data and corresponding contour boundary.
figure (8)
hold on
box on
set(gcf,'DefaultAxesFontSize',12);
p = plot(T_R,Hs_R,'k');
plot(T, Hs, '.')
title('Extreme Sea State Contour')
xlabel(TPlotLabel) 
ylabel('Significant Wave Height   H_s   [m]');
for i = 1:length(Time_r)
    [multmax,indmultmax] = max(Hs_R(:,i).*T_R(:,i));
    Hsmax = Hs_R(indmultmax,i);
    T_Hsmax = T_R(indmultmax,i);
    text(T_Hsmax,Hsmax,sprintf('%d yr',Time_r(i)),'FontWeight','bold','Rotation',-50,'BackgroundColor','w');
end
event_axis = axis;
hold off

%% If input is provided by user, include plots with steepness:
if exist('SteepMax','var')
    
% Figure 9: Plot steepness curve with extreme sea state contour
figure (9)
hold on
box on
set(gcf,'DefaultAxesFontSize',12);
plot(T_R,Hs_R,'k');
p = plot( T_vals(:,1), SteepH, 'r--', 'LineWidth', 1.2);
plot(T, Hs, '.')
title('Extreme Sea State Contour')
xlabel(TPlotLabel);
ylabel('Significant Wave Height   H_s   [m]');
for i = 1:length(Time_r)
    [multmax,indmultmax] = max(Hs_R(:,i).*T_R(:,i));
    Hsmax = Hs_R(indmultmax,i);
    T_Hsmax = T_R(indmultmax,i);
    text(T_Hsmax,Hsmax,sprintf('%d yr',Time_r(i)),'FontWeight','bold','Rotation',-50,'BackgroundColor','w');
end
axis(event_axis)
legend([p],'Wave breaking steepness');
hold off

% Figure 10: Plot data and contour with steepness included in contour.
figure (10)
hold on
box on
set(gcf,'DefaultAxesFontSize',12);
plot(T_R,Hs_R_2,'k');
p = plot( T_vals(:,1), SteepH, 'r--', 'LineWidth', 1.2);
plot(T, Hs, '.')
title('Extreme Sea State Contour - Including Steepness Limit')
xlabel(TPlotLabel);
ylabel('Significant Wave Height   H_s   [m]');
for i = 1:length(Time_r)
    [multmax,indmultmax] = max(Hs_R_2(:,i).*T_R(:,i));
    Hsmax = Hs_R_2(indmultmax,i);
    T_Hsmax = T_R(indmultmax,i);
    text(T_Hsmax,Hsmax,sprintf('%d yr',Time_r(i)),'FontWeight','bold','Rotation',-50,'BackgroundColor','w');
end
axis(event_axis)
legend([p],'Wave breaking steepness');
hold off

else
    figure(9); figure(10); % creates empty figures for saving later
end
%%
% Density Plots
% Figure 11: Plot of data density with extreme sea state contour. 
figure (11)
set(gcf,'DefaultAxesFontSize',12);
c = plot(T_R,Hs_R,'k');
hold on
box on
scatter(T,Hs,[],Count_val,'.')
xlabel(TPlotLabel);
ylabel('Significant Wave Height   H_s   [m]');
for i = 1:length(Time_r)
    [multmax,indmultmax] = max(Hs_R(:,i).*T_R(:,i));
    Hsmax = Hs_R(indmultmax,i);
    T_Hsmax = T_R(indmultmax,i);
    text(T_Hsmax,Hsmax,sprintf('%d yr',Time_r(i)),'FontWeight','bold','Rotation',-50,'BackgroundColor','w');
end
axis(event_axis)
t = colorbar('peer',gca);
set(get(t,'ylabel'),'string','Number of Points in Neighborhood')
hold off

% Figure 12: Plot of data density.
figure (12)
set(gcf,'DefaultAxesFontSize',12);
scatter(T,Hs,[],Count_val,'.')
box on
xlabel(TPlotLabel);
ylabel('Significant Wave Height   H_s   [m]');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'string','Number of Points in Neighborhood')

% If input is provided by user, include plot with steepness:

% Figure 13: Plot of data density with extreme sea state contour and 
% steepness curve. 
figure (13)
if exist('SteepMax','var')
set(gcf,'DefaultAxesFontSize',12);
c = plot(T_R,Hs_R,'k');
hold on
box on
p = plot( T_vals(:,1), SteepH, 'r--', 'LineWidth', 1.2);
legend([p],'Wave breaking steepness','Location','best');
scatter(T,Hs,[],Count_val,'.')
xlabel(TPlotLabel);
ylabel('Significant Wave Height   H_s   [m]');
for i = 1:length(Time_r)
    [multmax,indmultmax] = max(Hs_R_2(:,i).*T_R(:,i));
    Hsmax = Hs_R_2(indmultmax,i);
    T_Hsmax = T_R(indmultmax,i);
    text(T_Hsmax,Hsmax,sprintf('%d yr',Time_r(i)),'FontWeight','bold','Rotation',-50,'BackgroundColor','w');
end
axis(event_axis)
t = colorbar('peer',gca);
set(get(t,'ylabel'),'string','Number of Points in Neighborhood')
hold off
end

% Copyright 2015 Sandia Corporation. Under the terms of 
% Contract DE-AC04-94AL85000, there is a non-exclusive license for use of 
% this work by or on behalf of the U.S. Government. Export of this program 
% may require a license from the United States Government.