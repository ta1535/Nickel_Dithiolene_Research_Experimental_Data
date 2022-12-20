function Micro_GC_990_Data_Instant_Grapher_254_v02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify Path Address containing ASCII files of data
path_address = 'C:\Users\Halcy\Dropbox\Nickel_Dithiolene_Research_Experimental_Data\Raw_Data_Analysis\Experimental_Result\Data_14_Experiment\Compiled_Data\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change this number to look at specific GC Data file (counted from 1st
% file in directory, e.g. "5" would mean 5th GC graph found in a folder of
% multiple GC graph datas)
file_count = 88; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Searches for ASCII type files and file names
files_search = dir(strcat(path_address,'*.asc')); % specify file type here, e.g. *.asc is for ASCII type files
files_list = {files_search.name}; % Generates a string array list of all raw GC data graph files
file_directory_path_address = strings([1,length(files_list)]); % Generate list of file names for all GC graph data files
file_directory_path_address(file_count) = strcat(path_address,files_list(file_count)); % Generate string of full path address pointing to a specific GC data file
delimiterIn = ' '; % space separated between data values
GC_data = importdata(file_directory_path_address(file_count), delimiterIn);
GC_data_array_num = GC_data.data; % numerical portion of GC data
GC_data_array_text = GC_data.textdata; % string/text portion of GC data
Number_Points = str2double(extractBetween(replace(GC_data_array_text(14),sprintf('\t')," "),"Points: "," ")); % total number of data points read from text description

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparation/adjustment to raw data for peak analysis

% smooth raw data to rid of micro-peaks and noises
Channel_1_mV = smoothdata((GC_data_array_num(1:Number_Points)*(1e-5)),'gaussian',20); % Each data point's actual value is 10^5 smaller
Channel_2_mV = smoothdata((GC_data_array_num(Number_Points+1:end)*(1e-5)),'gaussian',10); % Each data point's actual value is 10^5 smaller
% note: channel 1 and 2 data is on the same vertical column; 1st half is channel 1 and 2nd half is channel 2.

%Channel_2_mV(6500:7000) = smoothdata(Channel_2_mV(6500:7000),'gaussian',100);

[Minimum_Channel_1,Min_1] = min(Channel_1_mV); % Find minimum mV value on Channel 1 for normalization purpose
[Minimum_Channel_2,Min_2] = min(Channel_2_mV); % Find minimum mV value on Channel 1 for normalization purpose
Channel_1_mV = Channel_1_mV - Minimum_Channel_1; % Data normalized so that all data points are positive with minimum value of 0
Channel_2_mV = Channel_2_mV - Minimum_Channel_2; % Data normalized so that all data points are positive with minimum value of 0

% each data point on on the time axis has a 0.01 second interval
Chromatograph_Time_Vector = [1:Number_Points]/100; % seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[GC_Raw_Data_1, GC_Raw_Data_2] = GC_Raw_Data_Analyzer(Chromatograph_Time_Vector,Channel_1_mV,Channel_2_mV);
[H2,O2,N2,C3H6,C3H8,CH4,CO2,CO,C2H4,C2H6,H2O,C4H10] = GC_Discrete_Analyzer(GC_Raw_Data_1,GC_Raw_Data_2);

% Channel 1 Summary:
Hydrogen_label = strcat("Hydrogen:   ",num2str(H2(end))," mV-s");
Oxygen_label = strcat("Oxygen:   ",num2str(O2(end))," mV-s");
Nitrogen_label = strcat("Nitrogen:   ",num2str(N2(end))," mV-s");
Methane_label = strcat("Methane:   ",num2str(CH4(end))," mV-s");
Carbon_Monoxide_label = strcat("Carbon Monoxide:   ",num2str(CO(end))," mV-s");

% Channel 2 Summary:
Carbon_Dioxide_label = strcat("Carbon Dioxide:   ",num2str(CO2(end))," mV-s");
Acetylene_Ethylene_label = strcat("Acetylene/Ethylene:   ",num2str(C2H4(end))," mV-s");
Ethane_label = strcat("Ethane:   ",num2str(C2H6(end))," mV-s");
Water_label = strcat("Water Vapor:   ",num2str(H2O(end))," mV-s");
Propylene_label = strcat("Propylene:   ",num2str(C3H6(end))," mV-s");
Propane_label = strcat("Propane:   ",num2str(C3H8(end))," mV-s");
Butane_label = strcat("Butane:   ",num2str(C4H10(end))," mV-s");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);

% Channel 1 Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1=subplot(2,1,1);
plot(Chromatograph_Time_Vector,Channel_1_mV,'b','linewidth',2,'markersize',15);
hold on;

% Hydrogen
plot(H2(3:4),H2(5:6),'r','linewidth',2,'markersize',15);
plot(H2(3),H2(5),'r|','linewidth',2,'markersize',15);
plot(H2(4),H2(6),'r|','linewidth',2,'markersize',15);
text(H2(1),H2(2),"H_2",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Oxygen
plot(O2(3:4),O2(5:6),'m','linewidth',2,'markersize',15);
plot(O2(3),O2(5),'m|','linewidth',2,'markersize',15);
plot(O2(4),O2(6),'m|','linewidth',2,'markersize',15);
text(O2(1),O2(2),"O_2",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Nitrogen
plot(N2(3:4),N2(5:6),'g','linewidth',2,'markersize',15);
plot(N2(3),N2(5),'g|','linewidth',2,'markersize',15);
plot(N2(4),N2(6),'g|','linewidth',2,'markersize',15);
text(N2(1),N2(2),"N_2",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Methane
plot(CH4(3:4),CH4(5:6),'c','linewidth',2,'markersize',15);
plot(CH4(3),CH4(5),'c|','linewidth',2,'markersize',15);
plot(CH4(4),CH4(6),'c|','linewidth',2,'markersize',15);
text(CH4(1),CH4(2),"CH_4",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Carbon monoxide
plot(CO(3:4),CO(5:6),'y','linewidth',2,'markersize',15);
plot(CO(3),CO(5),'y|','linewidth',2,'markersize',15);
plot(CO(4),CO(6),'y|','linewidth',2,'markersize',15);
text(CO(1),CO(2),"CO",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Summary table for Channel 1
dim1 = [.7 .5 .3 .3];
str1 = {Hydrogen_label;Oxygen_label;Nitrogen_label;Methane_label;Carbon_Monoxide_label};
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
legend('Channel 1 - Mol. Sieve 5A');

% Graph feature Channel 1
ylabel('GC unit [mV]');
xlim([0 Chromatograph_Time_Vector(end)]);
ylim([-max(Channel_1_mV)*0.2 max(Channel_1_mV)*1.2]);
set(gca,'FontName','Arial','FontSize',20,'linewidth',2,'TickLength',[0.025 0.025]);
pbaspect([3 1 1]);

% Channel 2 Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2=subplot(2,1,2);
plot(Chromatograph_Time_Vector,Channel_2_mV,'g','linewidth',2,'markersize',15);
hold on;

% Carbon dioxide
plot(CO2(3:4),CO2(5:6),'r','linewidth',2,'markersize',15);
plot(CO2(3),CO2(5),'r|','linewidth',2,'markersize',15);
plot(CO2(4),CO2(6),'r|','linewidth',2,'markersize',15);
text(CO2(1),CO2(2),"CO_2",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Ethylene/Acetylene
plot(C2H4(3:4),C2H4(5:6),'b','linewidth',2,'markersize',15);
plot(C2H4(3),C2H4(5),'b|','linewidth',2,'markersize',15);
plot(C2H4(4),C2H4(6),'b|','linewidth',2,'markersize',15);
text(C2H4(1),C2H4(2),"C_2H_4",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Ethane
plot(C2H6(3:4),C2H6(5:6),'m','linewidth',2,'markersize',15);
plot(C2H6(3),C2H6(5),'m|','linewidth',2,'markersize',15);
plot(C2H6(4),C2H6(6),'m|','linewidth',2,'markersize',15);
text(C2H6(1),C2H6(2),"C_2H_6",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Water Vapor
plot(H2O(3:4),H2O(5:6),'c','linewidth',2,'markersize',15);
plot(H2O(3),H2O(5),'c|','linewidth',2,'markersize',15);
plot(H2O(4),H2O(6),'c|','linewidth',2,'markersize',15);
text(H2O(1),H2O(2),"H_2O",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Propylene
plot(C3H6(3:4),C3H6(5:6),'b','linewidth',2,'markersize',15);
plot(C3H6(3),C3H6(5),'b|','linewidth',2,'markersize',15);
plot(C3H6(4),C3H6(6),'b|','linewidth',2,'markersize',15);
text(C3H6(1),C3H6(2),"C_3H_6",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Propane
plot(C3H8(3:4),C3H8(5:6),'m','linewidth',2,'markersize',15);
plot(C3H8(3),C3H8(5),'m|','linewidth',2,'markersize',15);
plot(C3H8(4),C3H8(6),'m|','linewidth',2,'markersize',15);
text(C3H8(1),C3H8(2),"C_3H_8",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Butane
plot(C4H10(3:4),C4H10(5:6),'r','linewidth',2,'markersize',15);
plot(C4H10(3),C4H10(5),'r|','linewidth',2,'markersize',15);
plot(C4H10(4),C4H10(6),'r|','linewidth',2,'markersize',15);
text(C4H10(1),C4H10(2),"C_4H_1_0",'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);

% Summary table for Channel 2
dim2 = [.7 .1 .3 .3];
str2 = {Carbon_Dioxide_label;Acetylene_Ethylene_label;Ethane_label;Water_label;Propylene_label;Propane_label;Butane_label};
annotation('textbox',dim2,'String',str2,'FitBoxToText','on');
legend('Channel 2 - PoraPLOT Q');

% Graph feature Channel 2
ylabel('GC unit [mV]');
xlabel('Retention Time [s]')
xlim([0 Chromatograph_Time_Vector(end)]);
ylim([-max(Channel_2_mV)*0.2 max(Channel_2_mV)*1.2]);
set(gca,'FontName','Arial','FontSize',20,'linewidth',2,'TickLength',[0.025 0.025]);
pbaspect([3 1 1]);

set(gcf, 'Position', [100, 100, 1200, 600]);
set(h1, 'Units', 'normalized');
set(h1, 'Position', [0.1,0.52,0.9,0.4]);
set(h2, 'Units', 'normalized');
set(h2, 'Position', [0.1,0.1,0.9,0.4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function [H2,O2,N2,C3H6,C3H8,CH4,CO2,CO,C2H4,C2H6,H2O,C4H10] = GC_Discrete_Analyzer(GC_Raw_Data_1,GC_Raw_Data_2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retention Time Parameters Input
% For Channel 1 MS5A Column
Channel_1_H2_RT_Lower_Bound = 18; % seconds
Channel_1_H2_RT_Upper_Bound = 20; % seconds
Channel_1_N2_RT_Lower_Bound = 29; % seconds
Channel_1_N2_RT_Upper_Bound = 35; % seconds
Channel_1_O2_RT_Lower_Bound = 22; % seconds
Channel_1_O2_RT_Upper_Bound = 25; % seconds
Channel_1_CH4_RT_Lower_Bound = 44; % seconds
Channel_1_CH4_RT_Upper_Bound = 47; % seconds
Channel_1_CO_RT_Lower_Bound = 60; % seconds
Channel_1_CO_RT_Upper_Bound = 90; % seconds

% For Channel 2 PoraplotQ Column
% Note: Ethylene/Acetylene are at same peak
Channel_2_C3H6_RT_Lower_Bound = 62; % seconds
Channel_2_C3H6_RT_Upper_Bound = 68.5; % seconds
Channel_2_C3H8_RT_Lower_Bound = 69; % seconds
Channel_2_C3H8_RT_Upper_Bound = 75; % seconds
Channel_2_C2H4_RT_Lower_Bound = 30; % seconds
Channel_2_C2H4_RT_Upper_Bound = 31; % seconds
Channel_2_CO2_RT_Lower_Bound = 26; % seconds
Channel_2_CO2_RT_Upper_Bound = 28; % seconds
Channel_2_C2H6_RT_Lower_Bound = 32.5; % seconds
Channel_2_C2H6_RT_Upper_Bound = 33.5; % seconds
Channel_2_H2O_RT_Lower_Bound = 40; % seconds
Channel_2_H2O_RT_Upper_Bound = 43; % seconds
Channel_2_C4H10_RT_Lower_Bound = 230; % seconds
Channel_2_C4H10_RT_Upper_Bound = 270; % seconds

% Medium Concentration RT
medium_limit = 5; % mV-s
Channel_2_C3H6_RT_Lower_Bound_medium = 64; % seconds
Channel_2_C3H6_RT_Upper_Bound_medium = 68; % seconds
Channel_2_C3H8_RT_Lower_Bound_medium = 70; % seconds
Channel_2_C3H8_RT_Upper_Bound_medium = 74; % seconds

% High Concentration RT
high_limit = 30; % mV-s
Channel_2_C3H6_RT_Lower_Bound_high = 52; % seconds
Channel_2_C3H6_RT_Upper_Bound_high = 54; % seconds
Channel_2_C3H8_RT_Lower_Bound_high = 56; % seconds
Channel_2_C3H8_RT_Upper_Bound_high = 58; % seconds

% Define variable vectors
Number_features = 7;
H2 = NaN(Number_features,1);
O2 = NaN(Number_features,1);
N2 = NaN(Number_features,1);
C3H6 = NaN(Number_features,1);
C3H8 = NaN(Number_features,1);
CH4 = NaN(Number_features,1);
CO2 = NaN(Number_features,1);
CO = NaN(Number_features,1);
C2H4 = NaN(Number_features,1);
C2H6 = NaN(Number_features,1);
H2O = NaN(Number_features,1);
C4H10 = NaN(Number_features,1);

H2(end) = 0;
O2(end) = 0;
N2(end) = 0;
C3H6(end) = 0;
C3H8(end) = 0;
CH4(end) = 0;
CO2(end) = 0;
CO(end) = 0;
C2H4(end) = 0;
C2H6(end) = 0;
H2O(end) = 0;
C4H10(end) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT vectors from raw datas
GC_Raw_Data_Channel_1_RT = GC_Raw_Data_1(:,1);
GC_Raw_Data_Channel_2_RT = GC_Raw_Data_2(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Label appropriate gas identity to peaks

count = 1;
while count <= length(GC_Raw_Data_Channel_1_RT)
    % Hydrogen check
    if GC_Raw_Data_Channel_1_RT(count) >= Channel_1_H2_RT_Lower_Bound && GC_Raw_Data_Channel_1_RT(count) <= Channel_1_H2_RT_Upper_Bound
        H2 = GC_Raw_Data_1(count,:);
    else % Oxygen check
        if GC_Raw_Data_Channel_1_RT(count) >= Channel_1_O2_RT_Lower_Bound && GC_Raw_Data_Channel_1_RT(count) <= Channel_1_O2_RT_Upper_Bound
            O2 = GC_Raw_Data_1(count,:);
        else % Nitrogen check
            if GC_Raw_Data_Channel_1_RT(count) >= Channel_1_N2_RT_Lower_Bound && GC_Raw_Data_Channel_1_RT(count) <= Channel_1_N2_RT_Upper_Bound
                N2 = GC_Raw_Data_1(count,:);
            else % Methane check
                if GC_Raw_Data_Channel_1_RT(count) >= Channel_1_CH4_RT_Lower_Bound && GC_Raw_Data_Channel_1_RT(count) <= Channel_1_CH4_RT_Upper_Bound
                    CH4 = GC_Raw_Data_1(count,:);
                else % Carbon monoxide check
                    if GC_Raw_Data_Channel_1_RT(count) >= Channel_1_CO_RT_Lower_Bound && GC_Raw_Data_Channel_1_RT(count) <= Channel_1_CO_RT_Upper_Bound
                        CO = GC_Raw_Data_1(count,:);
                    end
                end
            end
        end
    end

    count = count + 1;
end

count = 1;
propylene_count = 1;
while count <= length(GC_Raw_Data_Channel_2_RT)
    % Carbon dioxide check
    if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_CO2_RT_Lower_Bound && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_CO2_RT_Upper_Bound
        CO2 = GC_Raw_Data_2(count,:);
    else % Ethylene/Acetylene check
        if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_C2H4_RT_Lower_Bound && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_C2H4_RT_Upper_Bound
            C2H4 = GC_Raw_Data_2(count,:);
        else % Ethane check
            if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_C2H6_RT_Lower_Bound && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_C2H6_RT_Upper_Bound
                C2H6 = GC_Raw_Data_2(count,:);
            else % Propylene check
                if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_C3H6_RT_Lower_Bound && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_C3H6_RT_Upper_Bound
                    if propylene_count == 1
                        C3H6 = GC_Raw_Data_2(count,:);
                    end
                    propylene_count = propylene_count + 1;
                else % Propane check
                    if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_C3H8_RT_Lower_Bound && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_C3H8_RT_Upper_Bound
                        C3H8 = GC_Raw_Data_2(count,:);
                    else % Water check
                        if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_H2O_RT_Lower_Bound && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_H2O_RT_Upper_Bound
                            H2O = GC_Raw_Data_2(count,:);
                        else % Butane check
                            if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_C4H10_RT_Lower_Bound && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_C4H10_RT_Upper_Bound
                                C4H10 = GC_Raw_Data_2(count,:);
                            end
                        end
                    end
                end
            end
        end
    end

    % medium concentration propylene/propane
    if GC_Raw_Data_2(count,end) >= medium_limit && GC_Raw_Data_2(count,end) < high_limit
        if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_C3H6_RT_Lower_Bound_medium && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_C3H6_RT_Upper_Bound_medium
            C3H6 = GC_Raw_Data_2(count,:);
        else % Propane check
            if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_C3H8_RT_Lower_Bound_medium && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_C3H8_RT_Upper_Bound_medium
                C3H8 = GC_Raw_Data_2(count,:);
            end
        end
    end

    % high concentration propylene/propane
    if GC_Raw_Data_2(count,end) >= high_limit
        if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_C3H6_RT_Lower_Bound_high && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_C3H6_RT_Upper_Bound_high
            C3H6 = GC_Raw_Data_2(count,:);
        else % Propane check
            if GC_Raw_Data_Channel_2_RT(count) >= Channel_2_C3H8_RT_Lower_Bound_high && GC_Raw_Data_Channel_2_RT(count) <= Channel_2_C3H8_RT_Upper_Bound_high
                C3H8 = GC_Raw_Data_2(count,:);
            end
        end
    end

    count = count + 1;
end


end

function [GC_Raw_Data_1, GC_Raw_Data_2] = GC_Raw_Data_Analyzer(Chromatograph_Time_Vector,Channel_1_mV,Channel_2_mV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate slope, second derivative and instantaneous integral area (trapezoid rule) for each GC channel
% For Channel 1:
Channel_1_slope_count = 1;
while Channel_1_slope_count <= length(Channel_1_mV)
    if Channel_1_slope_count == 1
        dmV_dt_1 = 0;
        int_dt_1 = 0;
    else
        dmV_dt_1 = [dmV_dt_1; (Channel_1_mV(Channel_1_slope_count)-Channel_1_mV(Channel_1_slope_count-1))/0.01];
        int_dt_1 = [int_dt_1; (Channel_1_mV(Channel_1_slope_count)+Channel_1_mV(Channel_1_slope_count-1))/2*0.01];
    end
    Channel_1_slope_count = Channel_1_slope_count + 1;
end

% For Channel 2:
Channel_2_slope_count = 1;
while Channel_2_slope_count <= length(Channel_2_mV)
    if Channel_2_slope_count == 1
        dmV_dt_2 = 0;
        int_dt_2 = 0;
    else
        dmV_dt_2 = [dmV_dt_2; (Channel_2_mV(Channel_2_slope_count)-Channel_2_mV(Channel_2_slope_count-1))/0.01];
        int_dt_2 = [int_dt_2; (Channel_2_mV(Channel_2_slope_count)+Channel_2_mV(Channel_2_slope_count-1))/2*0.01];
    end
    Channel_2_slope_count = Channel_2_slope_count + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify peak locations on given GC data
[Channel_1_local_maxima_GC,Channel_1_maxima_index] = findpeaks(Channel_1_mV,'MinPeakWidth',0.2,'MinPeakProminence',0.01);
[Channel_2_local_maxima_GC,Channel_2_maxima_index] = findpeaks(Channel_2_mV,'MinPeakWidth',0.4,'MinPeakProminence',0.001);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Peak boundaries and peak area for each peak detected
% For Channel 1:
Channel_1_count = 1;
Channel_1_RT = NaN(1,length(Channel_1_maxima_index)); % keeps track of retention time for each peak maxima on channel 1
Channel_1_PA = NaN(1,length(Channel_1_maxima_index)); % keeps track of integrated peak area for each identified peak on channel 1
Channel_1_PA_left_index = NaN(1,length(Channel_1_maxima_index)); % keeps track of lower peak boundary for all detected peaks on channel 1
Channel_1_PA_right_index = NaN(1,length(Channel_1_maxima_index)); % keeps track of upper peak boundary for all detected peaks on channel 1
    
while Channel_1_count <= length(Channel_1_maxima_index)

    jump = 1;
    max_loop_limit = 500; % limit for max peak width (hundredths of second) per peak

    Channel_1_RT(Channel_1_count) = Chromatograph_Time_Vector(Channel_1_maxima_index(Channel_1_count)); 
    Channel_1_PA(Channel_1_count) = 0;
    Channel_1_PA_left_index(Channel_1_count) = Channel_1_maxima_index(Channel_1_count)-1;
    Channel_1_PA_right_index(Channel_1_count) = Channel_1_maxima_index(Channel_1_count)+1;

    % Identify peak boundaries using first derivative
    while Channel_1_PA_left_index(Channel_1_count) > (Channel_1_RT(Channel_1_count) - (max_loop_limit/100/2)) && ...
            dmV_dt_1(Channel_1_PA_left_index(Channel_1_count)) >= 1e-10 && ...
            Channel_1_PA_left_index(Channel_1_count) >= 1
        
        Channel_1_PA_left_index(Channel_1_count) = Channel_1_PA_left_index(Channel_1_count) - jump;
    end

    while Channel_1_PA_right_index(Channel_1_count) < (Channel_1_RT(Channel_1_count) + (max_loop_limit/100/2)) && ...
            dmV_dt_1(Channel_1_PA_right_index(Channel_1_count)) <= -1e-10 && ...
            Channel_1_PA_right_index(Channel_1_count) <= length(Chromatograph_Time_Vector)-1
        
        Channel_1_PA_right_index(Channel_1_count) = Channel_1_PA_right_index(Channel_1_count) + jump;
    end
    
    trigger_count = 10;
    
    % Fine-tune peak boundary points by keeping track of slope of secant line
    while trigger_count < 10

        % Shift peak boundary right
        Old_Slope = abs((Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))-Channel_1_mV(Channel_1_PA_left_index(Channel_1_count)))/...
            (Channel_1_PA_right_index(Channel_1_count)-Channel_1_PA_left_index(Channel_1_count)));

        Channel_1_PA_right_index(Channel_1_count) = Channel_1_PA_right_index(Channel_1_count) + jump;
        
        New_Slope = abs((Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))-Channel_1_mV(Channel_1_PA_left_index(Channel_1_count)))/...
            (Channel_1_PA_right_index(Channel_1_count)-Channel_1_PA_left_index(Channel_1_count)));

        while Old_Slope > New_Slope && Channel_1_PA_right_index(Channel_1_count) <= length(Chromatograph_Time_Vector)-1
            Channel_1_PA_right_index(Channel_1_count) = Channel_1_PA_right_index(Channel_1_count) + jump;
            Old_Slope = New_Slope;
            New_Slope = abs((Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))-Channel_1_mV(Channel_1_PA_left_index(Channel_1_count)))/...
            (Channel_1_PA_right_index(Channel_1_count)-Channel_1_PA_left_index(Channel_1_count)));
        end

        Channel_1_PA_right_index(Channel_1_count) = Channel_1_PA_right_index(Channel_1_count) - jump;

        % Shift peak boundary left
        Old_Slope = abs((Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))-Channel_1_mV(Channel_1_PA_left_index(Channel_1_count)))/...
            (Channel_1_PA_right_index(Channel_1_count)-Channel_1_PA_left_index(Channel_1_count)));

        Channel_1_PA_left_index(Channel_1_count) = Channel_1_PA_left_index(Channel_1_count) - jump;
        
        New_Slope = abs((Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))-Channel_1_mV(Channel_1_PA_left_index(Channel_1_count)))/...
            (Channel_1_PA_right_index(Channel_1_count)-Channel_1_PA_left_index(Channel_1_count)));

        while Old_Slope > New_Slope && Channel_1_PA_left_index(Channel_1_count) >= 1
            Channel_1_PA_left_index(Channel_1_count) = Channel_1_PA_left_index(Channel_1_count) - jump;
            Old_Slope = New_Slope;
            New_Slope = abs((Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))-Channel_1_mV(Channel_1_PA_left_index(Channel_1_count)))/...
            (Channel_1_PA_right_index(Channel_1_count)-Channel_1_PA_left_index(Channel_1_count)));
        end

        Channel_1_PA_left_index(Channel_1_count) = Channel_1_PA_left_index(Channel_1_count) + jump;

        trigger_count = trigger_count + 1;
    end

    trigger_count = 1;

    % Fine-tune peak boundary points by keeping track of updated area
    while trigger_count < 10

        % Shift peak boundary right
        Old_Area = sum(int_dt_1(Channel_1_PA_left_index(Channel_1_count):Channel_1_PA_right_index(Channel_1_count))) ...
            - (Channel_1_mV(Channel_1_PA_left_index(Channel_1_count))+Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))) ...
            *(Chromatograph_Time_Vector(Channel_1_PA_right_index(Channel_1_count))-Chromatograph_Time_Vector(Channel_1_PA_left_index(Channel_1_count)))/2;
        
        Channel_1_PA_right_index(Channel_1_count) = Channel_1_PA_right_index(Channel_1_count) + jump;
        
        New_Area = sum(int_dt_1(Channel_1_PA_left_index(Channel_1_count):Channel_1_PA_right_index(Channel_1_count))) ...
            - (Channel_1_mV(Channel_1_PA_left_index(Channel_1_count))+Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))) ...
            *(Chromatograph_Time_Vector(Channel_1_PA_right_index(Channel_1_count))-Chromatograph_Time_Vector(Channel_1_PA_left_index(Channel_1_count)))/2;
       
        while Old_Area < New_Area && Channel_1_PA_right_index(Channel_1_count) <= length(Chromatograph_Time_Vector)-1
            Channel_1_PA_right_index(Channel_1_count) = Channel_1_PA_right_index(Channel_1_count) + jump;
            Old_Area = New_Area;
            New_Area = sum(int_dt_1(Channel_1_PA_left_index(Channel_1_count):Channel_1_PA_right_index(Channel_1_count))) ...
            - (Channel_1_mV(Channel_1_PA_left_index(Channel_1_count))+Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))) ...
            *(Chromatograph_Time_Vector(Channel_1_PA_right_index(Channel_1_count))-Chromatograph_Time_Vector(Channel_1_PA_left_index(Channel_1_count)))/2;
        end

        Channel_1_PA_right_index(Channel_1_count) = Channel_1_PA_right_index(Channel_1_count) - jump;

        % Shift peak boundary left
        Old_Area = sum(int_dt_1(Channel_1_PA_left_index(Channel_1_count):Channel_1_PA_right_index(Channel_1_count))) ...
            - (Channel_1_mV(Channel_1_PA_left_index(Channel_1_count))+Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))) ...
            *(Chromatograph_Time_Vector(Channel_1_PA_right_index(Channel_1_count))-Chromatograph_Time_Vector(Channel_1_PA_left_index(Channel_1_count)))/2;
        
        Channel_1_PA_left_index(Channel_1_count) = Channel_1_PA_left_index(Channel_1_count) - jump;
        
        New_Area = sum(int_dt_1(Channel_1_PA_left_index(Channel_1_count):Channel_1_PA_right_index(Channel_1_count))) ...
            - (Channel_1_mV(Channel_1_PA_left_index(Channel_1_count))+Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))) ...
            *(Chromatograph_Time_Vector(Channel_1_PA_right_index(Channel_1_count))-Chromatograph_Time_Vector(Channel_1_PA_left_index(Channel_1_count)))/2;
        
        while Old_Area < New_Area && Channel_1_PA_left_index(Channel_1_count) >= 1
            Channel_1_PA_left_index(Channel_1_count) = Channel_1_PA_left_index(Channel_1_count) - jump;
            Old_Area= New_Area;
            New_Area = sum(int_dt_1(Channel_1_PA_left_index(Channel_1_count):Channel_1_PA_right_index(Channel_1_count))) ...
            - (Channel_1_mV(Channel_1_PA_left_index(Channel_1_count))+Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))) ...
            *(Chromatograph_Time_Vector(Channel_1_PA_right_index(Channel_1_count))-Chromatograph_Time_Vector(Channel_1_PA_left_index(Channel_1_count)))/2;
        end

        Channel_1_PA_left_index(Channel_1_count) = Channel_1_PA_left_index(Channel_1_count) + jump;

        trigger_count = trigger_count + 1;
    end

    Channel_1_PA(Channel_1_count) = sum(int_dt_1(Channel_1_PA_left_index(Channel_1_count):Channel_1_PA_right_index(Channel_1_count))) ...
            - (Channel_1_mV(Channel_1_PA_left_index(Channel_1_count))+Channel_1_mV(Channel_1_PA_right_index(Channel_1_count))) ...
            *(Chromatograph_Time_Vector(Channel_1_PA_right_index(Channel_1_count))-Chromatograph_Time_Vector(Channel_1_PA_left_index(Channel_1_count)))/2;

    Channel_1_count = Channel_1_count + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Peak boundaries and peak area for each peak detected
% For Channel 2:
high_area_limit = 50;
Channel_2_count = 1;
Channel_2_RT = NaN(1,length(Channel_2_maxima_index)); % keeps track of retention time for each peak maxima on channel 2
Channel_2_PA = NaN(1,length(Channel_2_maxima_index)); % keeps track of integrated peak area for each identified peak on channel 2
Channel_2_PA_left_index = NaN(1,length(Channel_2_maxima_index)); % keeps track of lower peak boundary for all detected peaks on channel 2
Channel_2_PA_right_index = NaN(1,length(Channel_2_maxima_index)); % keeps track of upper peak boundary for all detected peaks on channel 2
    
while Channel_2_count <= length(Channel_2_maxima_index)

    max_loop_limit = 500; % limit for max peak width (hundredths of second) per peak

    Channel_2_RT(Channel_2_count) = Chromatograph_Time_Vector(Channel_2_maxima_index(Channel_2_count)); 
    Channel_2_PA(Channel_2_count) = 0;
    Channel_2_PA_left_index(Channel_2_count) = Channel_2_maxima_index(Channel_2_count)-1;
    Channel_2_PA_right_index(Channel_2_count) = Channel_2_maxima_index(Channel_2_count)+1;

    % Identify peak boundaries using first derivative
    while Channel_2_PA_left_index(Channel_2_count) > (Channel_2_RT(Channel_2_count) - (max_loop_limit/100/2)) && ...
            dmV_dt_2(Channel_2_PA_left_index(Channel_2_count)) >= 1e-10 && ...
            Channel_2_PA_left_index(Channel_2_count) >= 1
        
        Channel_2_PA_left_index(Channel_2_count) = Channel_2_PA_left_index(Channel_2_count) - jump;
    end

    while Channel_2_PA_right_index(Channel_2_count) < (Channel_2_RT(Channel_2_count) + (max_loop_limit/100/2)) && ...
            dmV_dt_2(Channel_2_PA_right_index(Channel_2_count)) <= -1e-10 && ...
            Channel_2_PA_right_index(Channel_2_count) <= length(Chromatograph_Time_Vector)-1
        
        Channel_2_PA_right_index(Channel_2_count) = Channel_2_PA_right_index(Channel_2_count) + jump;
    end
    
    trigger_count = 1;

    % Fine-tune peak boundary points by keeping track of slope of secant line
    while trigger_count < 10

        % Shift peak boundary right
        Old_Slope = abs((Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))-Channel_2_mV(Channel_2_PA_left_index(Channel_2_count)))/...
            (Channel_2_PA_right_index(Channel_2_count)-Channel_2_PA_left_index(Channel_2_count)));

        Channel_2_PA_right_index(Channel_2_count) = Channel_2_PA_right_index(Channel_2_count) + jump;
        
        New_Slope = abs((Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))-Channel_2_mV(Channel_2_PA_left_index(Channel_2_count)))/...
            (Channel_2_PA_right_index(Channel_2_count)-Channel_2_PA_left_index(Channel_2_count)));

        while Old_Slope > New_Slope && Channel_2_PA_right_index(Channel_2_count) <= length(Chromatograph_Time_Vector)-1
            Channel_2_PA_right_index(Channel_2_count) = Channel_2_PA_right_index(Channel_2_count) + jump;
            Old_Slope = New_Slope;
            New_Slope = abs((Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))-Channel_2_mV(Channel_2_PA_left_index(Channel_2_count)))/...
            (Channel_2_PA_right_index(Channel_2_count)-Channel_2_PA_left_index(Channel_2_count)));
        end

        Channel_2_PA_right_index(Channel_2_count) = Channel_2_PA_right_index(Channel_2_count) - jump;

        % Shift peak boundary left
        Old_Slope = abs((Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))-Channel_2_mV(Channel_2_PA_left_index(Channel_2_count)))/...
            (Channel_2_PA_right_index(Channel_2_count)-Channel_2_PA_left_index(Channel_2_count)));

        Channel_2_PA_left_index(Channel_2_count) = Channel_2_PA_left_index(Channel_2_count) - jump;
        
        New_Slope = abs((Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))-Channel_2_mV(Channel_2_PA_left_index(Channel_2_count)))/...
            (Channel_2_PA_right_index(Channel_2_count)-Channel_2_PA_left_index(Channel_2_count)));

        while Old_Slope > New_Slope && Channel_2_PA_left_index(Channel_2_count) >= 1
            Channel_2_PA_left_index(Channel_2_count) = Channel_2_PA_left_index(Channel_2_count) - jump;
            Old_Slope = New_Slope;
            New_Slope = abs((Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))-Channel_2_mV(Channel_2_PA_left_index(Channel_2_count)))/...
            (Channel_2_PA_right_index(Channel_2_count)-Channel_2_PA_left_index(Channel_2_count)));
        end

        Channel_2_PA_left_index(Channel_2_count) = Channel_2_PA_left_index(Channel_2_count) + jump;

        trigger_count = trigger_count + 1;
    end

    trigger_count = 1;

    % Fine-tune peak boundary points by keeping track of updated area
    while trigger_count < 10

        % Shift peak boundary right
        Old_Area = sum(int_dt_2(Channel_2_PA_left_index(Channel_2_count):Channel_2_PA_right_index(Channel_2_count))) ...
            - (Channel_2_mV(Channel_2_PA_left_index(Channel_2_count))+Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))) ...
            *(Chromatograph_Time_Vector(Channel_2_PA_right_index(Channel_2_count))-Chromatograph_Time_Vector(Channel_2_PA_left_index(Channel_2_count)))/2;
        
        Channel_2_PA_right_index(Channel_2_count) = Channel_2_PA_right_index(Channel_2_count) + jump;
        
        New_Area = sum(int_dt_2(Channel_2_PA_left_index(Channel_2_count):Channel_2_PA_right_index(Channel_2_count))) ...
            - (Channel_2_mV(Channel_2_PA_left_index(Channel_2_count))+Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))) ...
            *(Chromatograph_Time_Vector(Channel_2_PA_right_index(Channel_2_count))-Chromatograph_Time_Vector(Channel_2_PA_left_index(Channel_2_count)))/2;
       
        while Old_Area < New_Area && Channel_2_PA_right_index(Channel_2_count) <= length(Chromatograph_Time_Vector)-1
            Channel_2_PA_right_index(Channel_2_count) = Channel_2_PA_right_index(Channel_2_count) + jump;
            Old_Area = New_Area;
            New_Area = sum(int_dt_2(Channel_2_PA_left_index(Channel_2_count):Channel_2_PA_right_index(Channel_2_count))) ...
            - (Channel_2_mV(Channel_2_PA_left_index(Channel_2_count))+Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))) ...
            *(Chromatograph_Time_Vector(Channel_2_PA_right_index(Channel_2_count))-Chromatograph_Time_Vector(Channel_2_PA_left_index(Channel_2_count)))/2;
        end

        Channel_2_PA_right_index(Channel_2_count) = Channel_2_PA_right_index(Channel_2_count) - jump;

        % Shift peak boundary left
        Old_Area = sum(int_dt_2(Channel_2_PA_left_index(Channel_2_count):Channel_2_PA_right_index(Channel_2_count))) ...
            - (Channel_2_mV(Channel_2_PA_left_index(Channel_2_count))+Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))) ...
            *(Chromatograph_Time_Vector(Channel_2_PA_right_index(Channel_2_count))-Chromatograph_Time_Vector(Channel_2_PA_left_index(Channel_2_count)))/2;
        
        Channel_2_PA_left_index(Channel_2_count) = Channel_2_PA_left_index(Channel_2_count) - jump;
        
        New_Area = sum(int_dt_2(Channel_2_PA_left_index(Channel_2_count):Channel_2_PA_right_index(Channel_2_count))) ...
            - (Channel_2_mV(Channel_2_PA_left_index(Channel_2_count))+Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))) ...
            *(Chromatograph_Time_Vector(Channel_2_PA_right_index(Channel_2_count))-Chromatograph_Time_Vector(Channel_2_PA_left_index(Channel_2_count)))/2;
        
        while Old_Area < New_Area && Channel_2_PA_left_index(Channel_2_count) >= 1
            Channel_2_PA_left_index(Channel_2_count) = Channel_2_PA_left_index(Channel_2_count) - jump;
            Old_Area= New_Area;
            New_Area = sum(int_dt_2(Channel_2_PA_left_index(Channel_2_count):Channel_2_PA_right_index(Channel_2_count))) ...
            - (Channel_2_mV(Channel_2_PA_left_index(Channel_2_count))+Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))) ...
            *(Chromatograph_Time_Vector(Channel_2_PA_right_index(Channel_2_count))-Chromatograph_Time_Vector(Channel_2_PA_left_index(Channel_2_count)))/2;
        end

        Channel_2_PA_left_index(Channel_2_count) = Channel_2_PA_left_index(Channel_2_count) + jump;

        trigger_count = trigger_count + 1;
    end

    Channel_2_PA(Channel_2_count) = sum(int_dt_2(Channel_2_PA_left_index(Channel_2_count):Channel_2_PA_right_index(Channel_2_count))) ...
            - (Channel_2_mV(Channel_2_PA_left_index(Channel_2_count))+Channel_2_mV(Channel_2_PA_right_index(Channel_2_count))) ...
            *(Chromatograph_Time_Vector(Channel_2_PA_right_index(Channel_2_count))-Chromatograph_Time_Vector(Channel_2_PA_left_index(Channel_2_count)))/2;

    if Channel_2_PA(Channel_2_count) > high_area_limit
        Channel_2_PA(Channel_2_count) = sum(int_dt_2(Channel_2_PA_left_index(Channel_2_count):Channel_2_PA_right_index(Channel_2_count)));
        Channel_2_mV(Channel_2_PA_left_index(Channel_2_count)) = 0;
        Channel_2_mV(Channel_2_PA_right_index(Channel_2_count)) = 0;
    end

    Channel_2_count = Channel_2_count + 1;
end

% 1. Retention time
% 2. Peak height
% 3. Left boundary retention time
% 4. Right boundary retention time
% 5. Left boundary mV
% 6. Right boundary mV
% 7. Peak area

GC_Raw_Data_1 = transpose([Channel_1_RT; transpose(Channel_1_local_maxima_GC); Chromatograph_Time_Vector(Channel_1_PA_left_index); Chromatograph_Time_Vector(Channel_1_PA_right_index);...
    transpose(Channel_1_mV(Channel_1_PA_left_index)); transpose(Channel_1_mV(Channel_1_PA_right_index)); Channel_1_PA]);

GC_Raw_Data_2 = transpose([Channel_2_RT; transpose(Channel_2_local_maxima_GC); Chromatograph_Time_Vector(Channel_2_PA_left_index); Chromatograph_Time_Vector(Channel_2_PA_right_index);...
    transpose(Channel_2_mV(Channel_2_PA_left_index)); transpose(Channel_2_mV(Channel_2_PA_right_index)); Channel_2_PA]);

end