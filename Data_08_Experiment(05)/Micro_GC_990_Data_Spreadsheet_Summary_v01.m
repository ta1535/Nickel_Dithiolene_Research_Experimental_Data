function Micro_GC_990_Data_Spreadsheet_Summary_v01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify Path Address containing ASCII files of data
path_address = 'C:\Users\Halcy\Dropbox\Nickel_Dithiolene_Research_Experimental_Data\Raw_Data_Analysis\Experimental_Result\Data_08_Experiment\Compiled_Data\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Searches for ASCII type files and file names
file_count = 1; 
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
% Output data for individual gas

Number_features = 7;

Hydrogen = NaN(Number_features,length(files_list));
Oxygen = NaN(Number_features,length(files_list));
Nitrogen = NaN(Number_features,length(files_list));
Propylene = NaN(Number_features,length(files_list));
Propane = NaN(Number_features,length(files_list));
Methane = NaN(Number_features,length(files_list));
Carbon_Dioxide = NaN(Number_features,length(files_list));
Carbon_Monoxide = NaN(Number_features,length(files_list));
Acetylene_Ethylene = NaN(Number_features,length(files_list));
Ethane = NaN(Number_features,length(files_list));
Water = NaN(Number_features,length(files_list));
Butane = NaN(Number_features,length(files_list));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop GC Analysis for every ASCII file data

while file_count <= length(files_list)
    file_directory_path_address(file_count) = strcat(path_address,files_list(file_count));
    GC_data = importdata(file_directory_path_address(file_count), delimiterIn);
    GC_data_array_num = GC_data.data; % numerical portion of GC data
    GC_data_array_text = GC_data.textdata; % string/text portion of GC data
    Number_Points = str2double(extractBetween(replace(GC_data_array_text(14),sprintf('\t')," "),"Points: "," ")); % total number of data points read from text description

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preparation/adjustment to raw data for peak analysis
    
    % smooth raw data to rid of micro-peaks and noises
    Channel_1_mV = smoothdata((GC_data_array_num(1:Number_Points)*(1e-5)),'gaussian',20); % Each data point's actual value is 10^5 smaller
    Channel_2_mV = smoothdata((GC_data_array_num(Number_Points+1:end)*(1e-5)),'gaussian',10); % Each data point's actual value is 10^5 smaller
    % note: channel 1 and 2 data is on the same vertical column; 1st half is channel 1 and 2nd half is channel 2.

    [Minimum_Channel_1,Min_1] = min(Channel_1_mV); % Find minimum mV value on Channel 1 for normalization purpose
    [Minimum_Channel_2,Min_2] = min(Channel_2_mV); % Find minimum mV value on Channel 1 for normalization purpose
    Channel_1_mV = Channel_1_mV - Minimum_Channel_1; % Data normalized so that all data points are positive with minimum value of 0
    Channel_2_mV = Channel_2_mV - Minimum_Channel_2; % Data normalized so that all data points are positive with minimum value of 0

    % each data point on on the time axis has a 0.01 second interval
    Chromatograph_Time_Vector = [1:Number_Points]/100; % seconds

    % Determine peak locations and area
    [GC_Raw_Data_1, GC_Raw_Data_2] = GC_Raw_Data_Analyzer(Chromatograph_Time_Vector,Channel_1_mV,Channel_2_mV);
    
    % Identify Individual Gas
    [H2,O2,N2,C3H6,C3H8,CH4,CO2,CO,C2H4,C2H6,H2O,C4H10] = GC_Discrete_Analyzer(GC_Raw_Data_1,GC_Raw_Data_2);
    
    Hydrogen(:,file_count) = H2;
    Oxygen(:,file_count) = O2;
    Nitrogen(:,file_count) = N2;
    Propylene(:,file_count) = C3H6;
    Propane(:,file_count) = C3H8;
    Methane(:,file_count) = CH4;
    Carbon_Dioxide(:,file_count) = CO2;
    Carbon_Monoxide(:,file_count) = CO;
    Acetylene_Ethylene(:,file_count) = C2H4;
    Ethane(:,file_count) = C2H6;
    Water(:,file_count) = H2O;
    Butane(:,file_count) = C4H10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine Time of Data Collection since initial collection

    time_stamp_year = str2double(extractBetween(files_list(file_count),1,4));
    time_stamp_month = str2double(extractBetween(files_list(file_count),6,7));
    time_stamp_day = str2double(extractBetween(files_list(file_count),9,10));
    time_stamp_hour = str2double(extractBetween(files_list(file_count),12,13));
    time_stamp_minute = str2double(extractBetween(files_list(file_count),15,16));
    time_stamp_second = str2double(extractBetween(files_list(file_count),18,19));
    
    time_stamp_current = datetime(time_stamp_year,time_stamp_month,time_stamp_day,time_stamp_hour,time_stamp_minute,time_stamp_second);
    
    if file_count == 1
        time_stamp_epoch = time_stamp_current;
    end

    data_collection_time_vector(file_count) = convertTo(time_stamp_current,'epochtime','Epoch',time_stamp_epoch);
    data_collection_time_vector = double(double(data_collection_time_vector));
    disp(strcat(num2str(file_count)," of ",num2str(length(files_list))," GC Data Files Processed."));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    file_count = file_count + 1;
    
    clear Channel_1_mV Channel_2_mV;
    clear Minimum_Channel_1 Minimum_Channel_2;
    clear GC_data GC_data_array_num GC_data_array_text;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare to export analyzed data to excel spreadsheet

RT_count = 1;
PH_count = 2;
LB_RT_count = 3;
RB_RT_count = 4;
LB_mV_count = 5;
RB_mV_count = 6;
PA_count = 7;

Data_Label = {'Time [seconds]', 'Hydrogen (Ch1)', 'Oxygen (Ch1)', 'Nitrogen (Ch1)', 'Methane (Ch1)', 'Carbon Monoxide (Ch1)', ...
    'Carbon Dioxide (Ch2)', 'Acetylene/Ethylene (Ch2)', 'Ethane (Ch2)', 'Water Vapor (Ch2)', 'Propylene (Ch2)', 'Propane (Ch2)', 'Butane (Ch2)'};

% Peak Area / Main data of interest (mV-s)
Summary_Data = [data_collection_time_vector; Hydrogen(PA_count,:); Oxygen(PA_count,:); Nitrogen(PA_count,:); Methane(PA_count,:); Carbon_Monoxide(PA_count,:); ...
    Carbon_Dioxide(PA_count,:); Acetylene_Ethylene(PA_count,:); Ethane(PA_count,:); Water(PA_count,:); Propylene(PA_count,:); Propane(PA_count,:); Butane(PA_count,:)];
Labeled_Table = array2table(Summary_Data.','VariableNames',Data_Label);

% Peak Retention Time (seconds)
RT = [data_collection_time_vector; Hydrogen(RT_count,:); Oxygen(RT_count,:); Nitrogen(RT_count,:); Methane(RT_count,:); Carbon_Monoxide(RT_count,:); ...
    Carbon_Dioxide(RT_count,:); Acetylene_Ethylene(RT_count,:); Ethane(RT_count,:); Water(RT_count,:); Propylene(RT_count,:); Propane(RT_count,:); Butane(RT_count,:)];
RT_Table = array2table(RT.','VariableNames',Data_Label);

% Peak Height (mV)
PH = [data_collection_time_vector; Hydrogen(PH_count,:); Oxygen(PH_count,:); Nitrogen(PH_count,:); Methane(PH_count,:); Carbon_Monoxide(PH_count,:); ...
    Carbon_Dioxide(PH_count,:); Acetylene_Ethylene(PH_count,:); Ethane(PH_count,:); Water(PH_count,:); Propylene(PH_count,:); Propane(PH_count,:); Butane(PH_count,:)];
PH_Table = array2table(PH.','VariableNames',Data_Label);

% Left boundary retention time (seconds)
LB_RT = [data_collection_time_vector; Hydrogen(LB_RT_count,:); Oxygen(LB_RT_count,:); Nitrogen(LB_RT_count,:); Methane(LB_RT_count,:); Carbon_Monoxide(LB_RT_count,:); ...
    Carbon_Dioxide(LB_RT_count,:); Acetylene_Ethylene(LB_RT_count,:); Ethane(LB_RT_count,:); Water(LB_RT_count,:); Propylene(LB_RT_count,:); Propane(LB_RT_count,:); Butane(LB_RT_count,:)];
LB_RT_Table = array2table(LB_RT.','VariableNames',Data_Label);

% Right boundary retention time (seconds)
RB_RT = [data_collection_time_vector; Hydrogen(RB_RT_count,:); Oxygen(RB_RT_count,:); Nitrogen(RB_RT_count,:); Methane(RB_RT_count,:); Carbon_Monoxide(RB_RT_count,:); ...
    Carbon_Dioxide(RB_RT_count,:); Acetylene_Ethylene(RB_RT_count,:); Ethane(RB_RT_count,:); Water(RB_RT_count,:); Propylene(RB_RT_count,:); Propane(RB_RT_count,:); Butane(RB_RT_count,:)];
RB_RT_Table = array2table(RB_RT.','VariableNames',Data_Label);

% Left boundary height (mV)
LB_mV = [data_collection_time_vector; Hydrogen(LB_mV_count,:); Oxygen(LB_mV_count,:); Nitrogen(LB_mV_count,:); Methane(LB_mV_count,:); Carbon_Monoxide(LB_mV_count,:); ...
    Carbon_Dioxide(LB_mV_count,:); Acetylene_Ethylene(LB_mV_count,:); Ethane(LB_mV_count,:); Water(LB_mV_count,:); Propylene(LB_mV_count,:); Propane(LB_mV_count,:); Butane(LB_mV_count,:)];
LB_mV_Table = array2table(LB_mV.','VariableNames',Data_Label);

% Right boundary height (mV)
RB_mV = [data_collection_time_vector; Hydrogen(RB_mV_count,:); Oxygen(RB_mV_count,:); Nitrogen(RB_mV_count,:); Methane(RB_mV_count,:); Carbon_Monoxide(RB_mV_count,:); ...
    Carbon_Dioxide(RB_mV_count,:); Acetylene_Ethylene(RB_mV_count,:); Ethane(RB_mV_count,:); Water(RB_mV_count,:); Propylene(RB_mV_count,:); Propane(RB_mV_count,:); Butane(RB_mV_count,:)];
RB_mV_Table = array2table(RB_mV.','VariableNames',Data_Label);

% Export of arrays
writetable(Labeled_Table,'GC_Summary_Data.xls','Sheet','Summary_Data');
writetable(RT_Table,'GC_Data.xls','Sheet','RT_Table');
writetable(PH_Table,'GC_Data.xls','Sheet','PH_Table');
writetable(LB_RT_Table,'GC_Data.xls','Sheet','LB_RT_Table');
writetable(RB_RT_Table,'GC_Data.xls','Sheet','RB_RT_Table');
writetable(LB_mV_Table,'GC_Data.xls','Sheet','LB_mV_Table');
writetable(RB_mV_Table,'GC_Data.xls','Sheet','RB_mV_Table');

disp('GC Data Processing Complete!');

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
                    C3H6 = GC_Raw_Data_2(count,:);
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
    
    trigger_count = 10;

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