function MEA_Experimental_Result_v01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify where to import data from
CA_Data = table2array(readtable('CA_001_C01.txt','NumHeaderLines',1));
GC_Summary_Data = table2array(readtable('GC_Summary_Data.xls','NumHeaderLines',1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input following parameters
number_cycle = 4; % total number of experimental cycles
cycle_time = 4*3600; % total duration of a single experimental cycle, seconds
oxidation_duration = 30*60; % total duration of oxidative phase, seconds
oxidation_start_time = 0; % time at which oxidation starts, seconds
reduction_duration = 15*60; % total duration of reductive phase, seconds
reduction_start_time = 3*3600; % time at which reduction starts, seconds

flow_rate_flush_N2 = 5; % sccm
flow_rate_propylene = 20; % sccm
flow_rate_propane = 0; % sccm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion Factor mV-s vs. PPM (by volume) from Calibration measurements
Propylene_Calibration_Conversion_slope = 4.14888715101403/(1e4); % mV-s per ppm
Propylene_Calibration_Conversion_b = 0;
Propane_Calibration_Conversion_slope = 5.82921089530086/(1e4); % mV-s per ppm
Propane_Calibration_Conversion_b = 0;

% Acceptably accurate
Nitrogen_Calibration_Conversion = 0.000414664; % mV-s per ppm
Carbon_monoxide_Calibration_Conversion = 0.000452229; % mV-s per ppm
Carbon_dioxide_Calibration_Conversion = 0.00253324; % mV-s per ppm
Acetylene_Ethylene_Calibration_Conversion = 0.002709107; % mV-s per ppm
Ethane_Calibration_Conversion = 0.001030689; % mV-s per ppm
Butane_Calibration_Conversion = 0.001185885; % mV-s per ppm

% Not very accurate
Water_Calibration_Conversion = 9.502232395/(1e4); % mV-s per ppm
Oxygen_Calibration_Conversion = 3.7908/(1e4); % mV-s per ppm
Hydrogen_Calibration_Conversion = 0.0326/(1e4); % mV-s per ppm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other parameters
Diameter_MEA = 1.5875; % cm
Area_MEA = pi()/4*Diameter_MEA^2; % cm2
Faraday_constant = 96485; % Faraday Constant, C/mol
Standard_Temperature = 0 + 273.15; % Kelvin
Standard_Pressure = 0.986923; % atm (from 1 bar);
R = 8.20573660809596e-5; % Universal Gas Constant, m3-atm/K-mol
MW_propylene = 42.08; % molar mass of propylene, g/mol
total_gas_molar_flow_rate = Standard_Pressure*(flow_rate_flush_N2*(1e-6))/R/Standard_Temperature/60; % molar flow rate, mol/sec
total_gas_molar_flow_rate_mix_feed = Standard_Pressure*((flow_rate_propylene+flow_rate_propane)*(1e-6))/R/Standard_Temperature/60; % molar flow rate, mol/sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converted Concentration Values

time = GC_Summary_Data(:,1);
hydrogen = GC_Summary_Data(:,2)/(Hydrogen_Calibration_Conversion);
oxygen = GC_Summary_Data(:,3)/(Oxygen_Calibration_Conversion);
nitrogen = GC_Summary_Data(:,4)/Nitrogen_Calibration_Conversion;
methane = GC_Summary_Data(:,5);
carbon_monoxide = GC_Summary_Data(:,6)/Carbon_monoxide_Calibration_Conversion;
carbon_dioxide = GC_Summary_Data(:,7)/Carbon_dioxide_Calibration_Conversion;
acetylene_ethylene = GC_Summary_Data(:,8)/Acetylene_Ethylene_Calibration_Conversion;
ethane = GC_Summary_Data(:,9)/Ethane_Calibration_Conversion;
water = GC_Summary_Data(:,10)/Water_Calibration_Conversion;
propylene = (GC_Summary_Data(:,11)-Propylene_Calibration_Conversion_b)/(Propylene_Calibration_Conversion_slope);
propane = (GC_Summary_Data(:,12)-Propane_Calibration_Conversion_b)/(Propane_Calibration_Conversion_slope);
butane = GC_Summary_Data(:,13)/Butane_Calibration_Conversion;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze and quantify separated propylene from GC Data

[reference_all,fit_natural,fit_combined_ECM,fit_deconvoluted_ECM,released_propylene_area] = ...
    Propylene_Separation_Analyzer(time,propylene,number_cycle,cycle_time,reduction_duration,reduction_start_time);

released_propylene_mole = released_propylene_area*total_gas_molar_flow_rate/(1e6);
released_propylene_mass = released_propylene_area*total_gas_molar_flow_rate*MW_propylene/1000/(1e6); % mass of purified propylene, kg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze and quantify electrochemical charges

[Oxidation_plot_points,Reductive_plot_points,Neutral_plot_points,Oxidation_Summary,Reduction_Summary] = ...
    CA_Analyzer(CA_Data,number_cycle,cycle_time,oxidation_duration,oxidation_start_time,reduction_duration,reduction_start_time,released_propylene_mass);

% During Reductive Phase:
time_start_H2 = reduction_start_time;
time_end_H2 = reduction_start_time + reduction_duration;
H2_area = simple_gas_area_integration(time,hydrogen,number_cycle,cycle_time,time_start_H2,time_end_H2);
released_H2_mole = H2_area*total_gas_molar_flow_rate/(1e6);

time_start_CO_reduction = reduction_start_time;
time_end_CO_reduction = reduction_start_time + reduction_duration;
CO2_area_reduction = simple_gas_area_integration(time,carbon_dioxide,number_cycle,cycle_time,time_start_CO_reduction,time_end_CO_reduction);
released_CO2_mole_reduction = CO2_area_reduction*total_gas_molar_flow_rate/(1e6);

% During Oxidative Phase:
time_start_CO2 = oxidation_start_time;
time_end_CO2 = oxidation_start_time + oxidation_duration;
CO2_area_oxidation = simple_gas_area_integration(time,carbon_dioxide,number_cycle,cycle_time,time_start_CO2,time_end_CO2);
released_CO2_mole_oxidation = CO2_area_oxidation*total_gas_molar_flow_rate_mix_feed/(1e6);

time_start_CO_oxidation = oxidation_start_time;
time_end_CO_oxidation = oxidation_start_time + oxidation_duration;
CO_area_oxidation = simple_gas_area_integration(time,carbon_monoxide,number_cycle,cycle_time,time_start_CO_oxidation,time_end_CO_oxidation);
released_CO_mole_oxidation = CO_area_oxidation*total_gas_molar_flow_rate_mix_feed/(1e6);

time_start_H2O = oxidation_start_time;
time_end_H2O = oxidation_start_time + oxidation_duration;
time_start_ref_H2O = 90*60;
H2O_area = water_vapor_area_integration(time,water,number_cycle,cycle_time,time_start_H2O,time_end_H2O,time_start_ref_H2O);
released_H2O_mole = H2O_area*total_gas_molar_flow_rate_mix_feed/(1e6);

% Faradaic Efficiency
% Reductive Phase:
FE_propylene_reduction = released_propylene_mole*2./(Reduction_Summary(2,:)/Faraday_constant);
FE_H2_reduction = released_H2_mole*2./(Reduction_Summary(2,:)/Faraday_constant);
FE_CO2_reduction = released_CO2_mole_reduction./(Reduction_Summary(2,:)/Faraday_constant);

% Oxidative Phase:
FE_CO_oxidation = released_CO_mole_oxidation./(Oxidation_Summary(2,:)/Faraday_constant);
FE_CO2_oxidation = released_CO2_mole_oxidation*2./(Oxidation_Summary(2,:)/Faraday_constant);
FE_H2O_oxidation = released_H2O_mole*2./(Oxidation_Summary(2,:)/Faraday_constant);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph Result for Figure #1
figure(1);
time = time/60; % convert seconds to minutes
x1 = 0; % domain limit of first graph subplot
x2 = number_cycle*cycle_time/60; % domain limit of first graph subplot
y1_propane_propylene = 0; % range of first graph subplot
y2_propane_propylene = 350; % range of first graph subplot
y1_water = 0; % range of second graph subplot
y2_water = 25000; % range of second graph subplot
y1_carbon = 0; % range of third graph subplot
y2_carbon = 100; % range of third graph subplot
CA_y_1 = -40; % range of fourth graph subplot
CA_y_2 = 40; % range of fourth graph subplot

% Subplot 1-1: Propylene and Propane
h1 = subplot(4,1,1);
plot(time,propylene,'.b','linewidth',2,'markersize',15);
hold on;
plot(time,propane,'.r','linewidth',2,'markersize',15);
graph_A_count = 1;
while graph_A_count <= number_cycle
    hold on;
    a_ox = area([((graph_A_count-1)*cycle_time+oxidation_start_time)/60 ((graph_A_count-1)*cycle_time+oxidation_start_time+oxidation_duration)/60],[max(propylene); max(propylene)]);
    a_ox.FaceColor = 'r';
    a_ox.FaceAlpha = 0.1;
    a_re = area([((graph_A_count-1)*cycle_time+reduction_start_time)/60 ((graph_A_count-1)*cycle_time+reduction_start_time+reduction_duration)/60],[max(propylene); max(propylene)]);
    a_re.FaceColor = 'b';
    a_re.FaceAlpha = 0.1;
    graph_A_count = graph_A_count + 1;
end
xlim([x1 x2]);
ylim([y1_propane_propylene y2_propane_propylene]);
legend('Propylene','Propane','Oxidative Phase','Reductive Phase','location','northeastoutside');
ylabel('Concentration [ppm]');
set(gca,'FontName','Arial','FontSize',10,'linewidth',2,'TickLength',[0.025 0.025],'xticklabel',[]);
pbaspect([8 1 1]);

% Subplot 1-2: Hydrogen, Oxygen, and Water
h2 = subplot(4,1,2);
plot(time,hydrogen,'.b','linewidth',2,'markersize',15);
hold on;
plot(time,oxygen,'.r','linewidth',2,'markersize',15);
plot(time,water,'.c','linewidth',2,'markersize',15);
graph_A_count = 1;
while graph_A_count <= number_cycle
    hold on;
    a_ox = area([((graph_A_count-1)*cycle_time+oxidation_start_time)/60 ((graph_A_count-1)*cycle_time+oxidation_start_time+oxidation_duration)/60],[1e7; 1e7]);
    a_ox.FaceColor = 'r';
    a_ox.FaceAlpha = 0.1;
    a_re = area([((graph_A_count-1)*cycle_time+reduction_start_time)/60 ((graph_A_count-1)*cycle_time+reduction_start_time+reduction_duration)/60],[1e7; 1e7]);
    a_re.FaceColor = 'b';
    a_re.FaceAlpha = 0.1;
    graph_A_count = graph_A_count + 1;
end
xlim([x1 x2]);
ylim([y1_water y2_water]);
ylabel('Concentration [ppm]');
legend('Hydrogen','Oxygen','Water Vapor','location','northeastoutside');
set(gca,'FontName','Arial','FontSize',10,'linewidth',2,'TickLength',[0.025 0.025],'xticklabel',[]);
pbaspect([8 1 1]);

% Subplot 1-3: Carbon monoxide and Carbon dioxide
h3 = subplot(4,1,3);
plot(time,carbon_monoxide,'.k','linewidth',2,'markersize',15);
hold on;
plot(time,carbon_dioxide,'.m','linewidth',2,'markersize',15);
graph_A_count = 1;
while graph_A_count <= number_cycle
    hold on;
    a_ox = area([((graph_A_count-1)*cycle_time+oxidation_start_time)/60 ((graph_A_count-1)*cycle_time+oxidation_start_time+oxidation_duration)/60],[1e7; 1e7]);
    a_ox.FaceColor = 'r';
    a_ox.FaceAlpha = 0.1;
    a_re = area([((graph_A_count-1)*cycle_time+reduction_start_time)/60 ((graph_A_count-1)*cycle_time+reduction_start_time+reduction_duration)/60],[1e7; 1e7]);
    a_re.FaceColor = 'b';
    a_re.FaceAlpha = 0.1;
    graph_A_count = graph_A_count + 1;
end
xlim([x1 x2]);
ylim([y1_carbon y2_carbon]);
ylabel('Concentration [ppm]');
legend('Carbon Monoxide','Carbon Dioxide','location','northeastoutside');
set(gca,'FontName','Arial','FontSize',10,'linewidth',2,'TickLength',[0.025 0.025],'xticklabel',[]);
pbaspect([8 1 1]);

% Subplot 1-4: Chronoamperometry
h4 = subplot(4,1,4);
plot(Oxidation_plot_points(1,:)/60,Oxidation_plot_points(2,:)/Area_MEA,'r','linewidth',2,'markersize',15);
hold on;
plot(Reductive_plot_points(1,:)/60,Reductive_plot_points(2,:)/Area_MEA,'b','linewidth',2,'markersize',15);
plot(Neutral_plot_points(1,:)/60,Neutral_plot_points(2,:)/Area_MEA,'g','linewidth',2,'markersize',15);
graph_B_count = 1;
while graph_B_count <= number_cycle
    hold on;
    a_ox = area([((graph_B_count-1)*cycle_time+oxidation_start_time)/60 ((graph_B_count-1)*cycle_time+oxidation_start_time+oxidation_duration)/60],[max(propylene); max(propylene)]);
    a_ox.FaceColor = 'r';
    a_ox.FaceAlpha = 0.1;
    a_ox2 = area([((graph_B_count-1)*cycle_time+oxidation_start_time)/60 ((graph_B_count-1)*cycle_time+oxidation_start_time+oxidation_duration)/60],[-max(propylene); -max(propylene)]);
    a_ox2.FaceColor = 'r';
    a_ox2.FaceAlpha = 0.1;
    a_re = area([((graph_B_count-1)*cycle_time+reduction_start_time)/60 ((graph_B_count-1)*cycle_time+reduction_start_time+reduction_duration)/60],[max(propylene); max(propylene)]);
    a_re.FaceColor = 'b';
    a_re.FaceAlpha = 0.1;
    a_re2 = area([((graph_B_count-1)*cycle_time+reduction_start_time)/60 ((graph_B_count-1)*cycle_time+reduction_start_time+reduction_duration)/60],[-max(propylene); -max(propylene)]);
    a_re2.FaceColor = 'b';
    a_re2.FaceAlpha = 0.1;
    text(((graph_B_count-1)*cycle_time+oxidation_start_time)/60,CA_y_1/2,strcat('+',num2str(Oxidation_Summary(2,graph_B_count),'%9.3f'),' C'));
    text(((graph_B_count-1)*cycle_time+reduction_start_time)/60,CA_y_2/2,strcat('-',num2str(Reduction_Summary(2,graph_B_count),'%9.3f'),' C'));
    graph_B_count = graph_B_count + 1;
end
legend(strcat("+",num2str(mean(Oxidation_Summary(1,:))),"V Oxidatative Current"),strcat(num2str(mean(Reduction_Summary(1,:))),"V Reductive Current"),"0V Applied",'location','northeastoutside');
xlim([x1 x2]);
ylim([CA_y_1 CA_y_2]);
ylabel('Current Density [mA/cm^2]');
xlabel('Time Elapsed [minutes]');
set(gca,'FontName','Arial','FontSize',10,'linewidth',2,'TickLength',[0.025 0.025]);
pbaspect([8 1 1]);

% Spacing/Positioning of individual subplots
set(h1, 'Units', 'normalized');
set(h1, 'Position', [0.05,0.6,0.8,0.4]);
set(h2, 'Units', 'normalized');
set(h2, 'Position', [0.05,0.4,0.8,0.4]);
set(h3, 'Units', 'normalized');
set(h3, 'Position', [0.05,0.2,0.8,0.4])
set(h4, 'Units', 'normalized');
set(h4, 'Position', [0.05,0,0.8,0.4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph Result for Figure #2
figure(2);
x_PS_1 = -2000; % domain limit of first graph subplot
x_PS_2 = 2000; % domain limit of first graph subplot
y_combined_1 = 0; % range limit of first graph subplot
y_combined_2 = 250; % range limit of first graph subplot
y_decon_1 = 0; % range limit of second graph subplot
y_decon_2 = 100; % range limit of second graph subplot
y_FE_1 = 0; % range limit of third graph subplot
y_FE_2 = 1; % range limit of third graph subplot

cycle_count = 1;
while cycle_count <= number_cycle
    h1 = subplot(3,number_cycle,cycle_count);
    plot(reference_all(:,(cycle_count-1)*2+1),reference_all(:,cycle_count*2),'.b','linewidth',2,'markersize',15,'HandleVisibility','off');
    hold on;
    plot(fit_natural(:,1),fit_natural(:,cycle_count+1),'g','linewidth',2,'markersize',15);
    plot(fit_combined_ECM(:,1),fit_combined_ECM(:,cycle_count+1),'r','linewidth',2,'markersize',15);
    xlim([x_PS_1 x_PS_2]);
    ylim([y_combined_1 y_combined_2]);
    ylabel('C_3H_6 Concentration [ppm]');
    xlabel('Time relative to P.S. [s]');
    set(gca,'FontName','Arial','FontSize',10,'linewidth',2,'TickLength',[0.025 0.025]);
    pbaspect([1 1 1]);
    
    h2 = subplot(3,number_cycle,cycle_count+number_cycle);
    plot(fit_deconvoluted_ECM(:,1),fit_deconvoluted_ECM(:,cycle_count+1),'r','linewidth',2,'markersize',15);
    xlim([x_PS_1 x_PS_2]);
    ylim([y_decon_1 y_decon_2]);
    ylabel('C_3H_6 Concentration [ppm]');
    xlabel('Time relative to P.S. [s]');
    set(gca,'FontName','Arial','FontSize',10,'linewidth',2,'TickLength',[0.025 0.025]);
    pbaspect([1 1 1]);
    
    cycle_count = cycle_count + 1;
end

h3 = subplot(3,1,3);
X_bar = 1:number_cycle;
Y_bar = transpose([FE_CO_oxidation;FE_CO2_oxidation;FE_propylene_reduction;FE_H2_reduction]);
bar_graph = bar(X_bar,Y_bar);
gradient_number = linspace(0,0.8,2);
bar_graph(1).FaceColor = [1 gradient_number(2) 0];
bar_graph(2).FaceColor = [1 gradient_number(1) 0];
bar_graph(3).FaceColor = [0 gradient_number(2) 1];
bar_graph(4).FaceColor = [0 gradient_number(1) 1];

xlabel('Experimental Run Cycle');
ylabel('Faradic Efficiency');
legend('CO','CO_2','C_3H_6','H_2','NumColumns',2,'location','northeastoutside');
leg = legend('show');
title(leg, 'Oxid. FE      Red. FE');
count = 1;
while count <= 4
    xtips = bar_graph(count).XEndPoints; ytips = bar_graph(count).YEndPoints;
    cycle_count = 1;
    while cycle_count <= number_cycle
        text(xtips(cycle_count),ytips(cycle_count),num2str(Y_bar(cycle_count,count),'%0.4f'),'vert','bottom','horiz','left','Fontsize',10,'rotation',60);
        cycle_count = cycle_count + 1;
    end
    count = count + 1;
end
ylim([y_FE_1 y_FE_2]);
set(gca,'YScale','log');
set(gca,'FontName','Arial','FontSize',10,'linewidth',2,'TickLength',[0.025 0.025]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data_Label = {'Propylene separated per unit area [mol/cm2]', 'FEO (H2O)', 'FEO (CO)', 'FEO (CO2)', 'FER (C3H6)', 'FER (H2)', 'FER (CO2)'};

final_output = transpose([released_propylene_mole;FE_H2O_oxidation;FE_CO_oxidation;FE_CO2_oxidation;FE_propylene_reduction;FE_H2_reduction;FE_CO2_reduction]);
final_output_table = array2table(final_output,'VariableNames',Data_Label);
writetable(final_output_table,'Processed_Data.xls','Sheet','Experimental Result');
writematrix(reference_all,'Processed_Data.xls','Sheet','reference_points');
writematrix(fit_natural,'Processed_Data.xls','Sheet','fit_natural');
writematrix(fit_combined_ECM,'Processed_Data.xls','Sheet','fit_combined_ECM');

export_file_address = "C:\Users\Halcy\Dropbox\Nickel_Dithiolene_Research_Experimental_Data\Data_Trends_Analysis\Deconvoluted_ECM_15.xls";
export_array = [fit_deconvoluted_ECM(:,1) fit_deconvoluted_ECM(:,3)];
writematrix(export_array,export_file_address);

end

function integrated_gas_area = simple_gas_area_integration(time,gas,number_cycle,cycle_time,time_start,time_end)

integrated_gas_area = zeros(1,number_cycle);
cycle_count = 1;
count = 1;

while cycle_count <= number_cycle
    sort_count = 1;
    while count <= length(gas) && time(count) <= cycle_count*cycle_time
        if time(count) >= (cycle_count-1)*cycle_time+time_start && time(count) <= (cycle_count-1)*cycle_time+time_end
            if sort_count == 1
                custom_gas = [time(count) gas(count)];
            else
                custom_gas = [custom_gas; time(count) gas(count)];
            end
            sort_count = sort_count + 1;
        end
        count = count + 1;
    end
    integrated_gas_area(cycle_count) = trapz(custom_gas(:,1),custom_gas(:,2));
    clear custom_gas;
    cycle_count = cycle_count + 1;
end

end

function water_area = water_vapor_area_integration(time,water,number_cycle,cycle_time,time_start,time_end,time_start_ref)

water_area = zeros(1,number_cycle);

cycle_count = 1;
count = 1;

while cycle_count <= number_cycle
    sort_count = 1;
    ref_count = 1;
    while count <= length(water) && time(count) <= cycle_count*cycle_time
        if time(count) >= (cycle_count-1)*cycle_time+time_start && time(count) <= (cycle_count-1)*cycle_time+time_end
            if sort_count == 1
                custom_water = [time(count) water(count)];
            else
                custom_water = [custom_water; time(count) water(count)];
            end
            sort_count = sort_count + 1;
        end

        if time(count) >= (cycle_count-1)*cycle_time+time_start_ref && time(count) <= (cycle_count-1)*cycle_time+time_start_ref+(time_end-time_start)
            if ref_count == 1
                custom_water_ref = [time(count) water(count)];
            else
                custom_water_ref = [custom_water_ref; time(count) water(count)];
            end
            ref_count = ref_count + 1;
        end
        count = count + 1;
    end
    water_area(cycle_count) = trapz(custom_water(:,1),custom_water(:,2))-trapz(custom_water_ref(:,1),custom_water_ref(:,2));

    clear custom_gas;
    clear custom_water_ref;
    cycle_count = cycle_count + 1;
end

end

function [reference_points_all,simulated_points_natural_desorption,simulated_points_combined_desorption,simulated_points_ECM_deconvoluted_desorption,released_propylene_area]...
    = Propylene_Separation_Analyzer(time,propylene,number_cycle,cycle_time,reduction_duration,reduction_start_time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range parameters when fitting curves onto actual data points

natural_desorption_curve_fitting_range_before = 60*60; % range of natural desorption time before reduction starts, seconds
natural_desorption_curve_fitting_range_after = 60*60; % range of natural desorption time after reduction ends, seconds

ECM_desorption_curve_fitting_range_before = 150; % range of ECM desorption time before reduction starts, seconds
ECM_desorption_curve_fitting_range_after = 15*60; % range of ECM desorption time after reduction ends, seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for generated curve fitting simulation (for graphing purpose)

number_points_fitting_curve = 1000;
start_time_sim = -1*natural_desorption_curve_fitting_range_before; % start time relative to P.S. for graphing display
end_time_sim = natural_desorption_curve_fitting_range_after; % end time relative to P.S. for graphing display
curve_time_vector = linspace(start_time_sim,end_time_sim,number_points_fitting_curve); % time vector axis to use for graphing display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort data between natural desorption points and points experiencing the
% effect of both electrochemical modulation (ECM) and natural desorption

% reference points are stored in pairs per cycle (e.g. 8 columns of data for 4 cycles)
reference_points_all = NaN(number_cycle*2,100);
released_propylene_area = zeros(1,number_cycle);

cycle_count = 1;
count = 1;

while cycle_count <= number_cycle

    natural_desorption_reference_points = [0 0]; % for natural desorption only
    natural_desorption_reference_count = 1;
    
    ECM_induced_desorption_reference_points = [0 0]; % for ECM combined with natural desorption
    ECM_induced_desorption_reference_count = 1;

    reference_count = 1;

    while count <= length(propylene) && time(count) <= cycle_count*cycle_time
        if time(count) >= cycle_time*(cycle_count-1)+reduction_start_time-natural_desorption_curve_fitting_range_before &&...
                time(count) <= cycle_time*(cycle_count-1)+reduction_start_time+reduction_duration+natural_desorption_curve_fitting_range_after
            reference_points_all((cycle_count-1)*2+1,reference_count) = time(count)-cycle_time*(cycle_count-1)-reduction_start_time;
            reference_points_all(cycle_count*2,reference_count) = propylene(count);
            reference_count = reference_count + 1;

            if time(count) >= cycle_time*(cycle_count-1)+reduction_start_time-ECM_desorption_curve_fitting_range_before &&...
                    time(count) <= cycle_time*(cycle_count-1)+reduction_start_time+reduction_duration+ECM_desorption_curve_fitting_range_after
                if ECM_induced_desorption_reference_count == 1
                    ECM_induced_desorption_reference_points(1) = time(count)-cycle_time*(cycle_count-1)-reduction_start_time;
                    ECM_induced_desorption_reference_points(2) = propylene(count);
                else
                    ECM_induced_desorption_reference_points = [ECM_induced_desorption_reference_points; (time(count)-cycle_time*(cycle_count-1)-reduction_start_time) propylene(count)];
                end
                ECM_induced_desorption_reference_count = ECM_induced_desorption_reference_count + 1;
            
            else
                if natural_desorption_reference_count == 1
                    natural_desorption_reference_points(1) = time(count)-cycle_time*(cycle_count-1)-reduction_start_time;
                    natural_desorption_reference_points(2) = propylene(count);
                else
                    natural_desorption_reference_points = [natural_desorption_reference_points; (time(count)-cycle_time*(cycle_count-1)-reduction_start_time) propylene(count)];
                end
                natural_desorption_reference_count = natural_desorption_reference_count + 1;
            
            end
        end
        count = count + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Polynomial curve fitting for natural desorption reference points

    propylene_natural_desorption_polyfit = polyfit(natural_desorption_reference_points(:,1),natural_desorption_reference_points(:,2),4);
    polyfitted_natural_desorption_y = polyval(propylene_natural_desorption_polyfit,curve_time_vector);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Increase number of points on rising slope to balance weighted fitting

    if ECM_induced_desorption_reference_count > 4

        check_count_max = 4;
        y_difference = zeros(1,check_count_max);
        check_count = 1;

        while check_count <= check_count_max
            y_difference(check_count) = ECM_induced_desorption_reference_points(check_count+1,2) - ECM_induced_desorption_reference_points(check_count,2);
            check_count = check_count + 1;
        end

        [max_diff,max_diff_index] = max(y_difference);

        add_points_propylene_x_1 = linspace(ECM_induced_desorption_reference_points(max_diff_index,1),ECM_induced_desorption_reference_points(max_diff_index+1,1),(length(ECM_induced_desorption_reference_points)-2));
        add_points_propylene_y_1 = linspace(ECM_induced_desorption_reference_points(max_diff_index,2),ECM_induced_desorption_reference_points(max_diff_index+1,2),(length(ECM_induced_desorption_reference_points)-2));
        add_points_propylene = [transpose(add_points_propylene_x_1) transpose(add_points_propylene_y_1)];
        
        ECM_induced_desorption_reference_points = [ECM_induced_desorption_reference_points; add_points_propylene];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deconvolute natural desorption from ECM reference points

    polyfitted_natural_desorption_y_remove = polyval(propylene_natural_desorption_polyfit,ECM_induced_desorption_reference_points(:,1));
    ECM_induced_desorption_reference_points_adjusted_y = ECM_induced_desorption_reference_points(:,2)-polyfitted_natural_desorption_y_remove;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nonlinear fitting of deconvoluted reference points with skewed Gaussian distribution curve 
    % b(1) = horizontal spread/distribution width of peak
    % b(2) = skewness of peak (e.g. positive for left leaning peak)
    % b(3) = height of peak
    % b(4) = peak location on x-axis

    model_function_ECM_desorption = @(b,x) b(3)./sqrt(2*pi()).*exp(-(x-b(4)).^2/(b(1)^2)).*(1+erf(b(2)*(x-b(4))./sqrt(2)));
    beta_0 = [1000 0.01 100 0];

    beta_ECM_desorption = nlinfit(ECM_induced_desorption_reference_points(:,1),ECM_induced_desorption_reference_points_adjusted_y,model_function_ECM_desorption,beta_0);
    fitted_deconvoluted_ECM_desorption_points_y = feval(model_function_ECM_desorption,beta_ECM_desorption,curve_time_vector);
    fitted_combined_desorption_points_y = fitted_deconvoluted_ECM_desorption_points_y + polyfitted_natural_desorption_y;

    natural_desorption_integrated_area = trapz(curve_time_vector,polyfitted_natural_desorption_y);
    combined_desorption_integrated_area = trapz(curve_time_vector,fitted_combined_desorption_points_y);
    released_propylene_area(cycle_count) = combined_desorption_integrated_area - natural_desorption_integrated_area;

    if cycle_count == 1
        simulated_points_natural_desorption = transpose([curve_time_vector;polyfitted_natural_desorption_y]);
        simulated_points_combined_desorption = transpose([curve_time_vector;fitted_combined_desorption_points_y]);
        simulated_points_ECM_deconvoluted_desorption = transpose([curve_time_vector;fitted_deconvoluted_ECM_desorption_points_y]);
    else
        % simulated natural desorption line contains time vector and
        % y-coordinates for each cycle (e.g. 5 columns of data for 4
        % cycles)
        simulated_points_natural_desorption = [simulated_points_natural_desorption, transpose(polyfitted_natural_desorption_y)]; 

        % simulated combined desorption line containes time vector and
        % y-coordinates for each cycle (e.g. 5 columns of data for 4
        % cycles)
        simulated_points_combined_desorption = [simulated_points_combined_desorption, transpose(fitted_combined_desorption_points_y)];

        % simulated deconvoluted ECM desorption line containes time vector
        % and y-coordinates for each cycle (e.g. 5 columns of data for 4
        % cycles)
        simulated_points_ECM_deconvoluted_desorption = [simulated_points_ECM_deconvoluted_desorption, transpose(fitted_deconvoluted_ECM_desorption_points_y)];
    end
    
    clear natural_desorption_reference_points;
    clear ECM_induced_desorption_reference_points;
    clear add_points_propylene_x_1;
    clear add_points_propylene_y_1;
    clear add_points_propylene;
    clear polyfitted_natural_desorption_y_remove;
    clear ECM_induced_desorption_reference_points_adjusted_y;
    clear fitted_deconvoluted_ECM_desorption_points_y;
    clear fitted_deconvoluted_ECM_desorption_points;
    clear fitted_combined_desorption_points_y;
    clear fitted_combined_desorption_points;
    clear natural_desorption_integrated_area;
    clear combined_desorption_integrated_area;

    cycle_count = cycle_count + 1;
end
reference_points_all = transpose(reference_points_all);
end

function [Oxidation_plot_points,Reductive_plot_points,Neutral_plot_points,Oxidation_Summary,Reduction_Summary] = ...
    CA_Analyzer(CA_Data,number_cycle,cycle_time,oxidation_duration,oxidation_start_time,reduction_duration,reduction_start_time,propylene_mass_separated)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepares Chronoamperometry data
CA_time = CA_Data(:,1);
CA_Potential = CA_Data(:,2);
CA_Current = CA_Data(:,3);
CA_Current_Smooth = smoothdata(CA_Current,'gaussian',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each feature stores different
% 1. Average Potential on given cycle (V)
% 2. Total Coloumbic Charge applied on given cycle (C)
% 3. Total Energy Applied on given cycle (J)
% 4. Total Energy per kg propylene separated (MJ/kg)

number_features = 4;
Oxidation_Summary = zeros(number_features,number_cycle);
Reduction_Summary = zeros(number_features,number_cycle);

Oxidation_plot_points = NaN(2,length(CA_Current_Smooth));
Reductive_plot_points = NaN(2,length(CA_Current_Smooth));
Neutral_plot_points = NaN(2,length(CA_Current_Smooth));

cycle_count = 1;
count = 1;

while cycle_count <= number_cycle
    oxidation_count = 1;
    reduction_count = 1;
    Oxidation_index_temp = [0 1];
    Reduction_index_temp = [0 1];
    Oxidation_Potential_temp = 0;
    Reduction_Potential_temp = 0;

    while count <= length(CA_time) && CA_time(count) <= cycle_count*cycle_time 
        if CA_Potential(count) > 0.01
            Oxidation_plot_points(1,count) = CA_time(count);
            Oxidation_plot_points(2,count) = CA_Current_Smooth(count);
            if oxidation_count == 1
                Oxidation_Potential_temp = CA_Potential(count);
                Oxidation_index_temp(1) = count;
            else
                Oxidation_Potential_temp = [Oxidation_Potential_temp CA_Potential(count)];
                Oxidation_index_temp(2) = count;
            end
            oxidation_count = oxidation_count + 1;
            reduction_count = 1;
        else
            if CA_Potential(count) < 0.01 && CA_Potential(count) > -0.01
                Neutral_plot_points(1,count) = CA_time(count);
                Neutral_plot_points(2,count) = CA_Current_Smooth(count);
            else
                if CA_Potential(count) < -0.01
                    Reductive_plot_points(1,count) = CA_time(count);
                    Reductive_plot_points(2,count) = CA_Current_Smooth(count);
                    if reduction_count == 1
                        Reduction_Potential_temp = CA_Potential(count);
                        Reduction_index_temp(1) = count;
                    else
                        Reduction_Potential_temp = [Reduction_Potential_temp CA_Potential(count)];
                        Reduction_index_temp(2) = count;
                    end
                    reduction_count = reduction_count + 1;
                    oxidation_count = 1;
                end
            end
        end
        count = count + 1;
    end

    Oxidation_Summary(1,cycle_count) = round(mean(Oxidation_Potential_temp),1); % Volt
    Oxidation_Summary(2,cycle_count) = trapz(CA_time(Oxidation_index_temp(1):Oxidation_index_temp(2)),CA_Current_Smooth(Oxidation_index_temp(1):Oxidation_index_temp(2)))/1000; % Coulomb
    Oxidation_Summary(3,cycle_count) = Oxidation_Summary(2,cycle_count)*Oxidation_Summary(1,cycle_count); % Joules
    Oxidation_Summary(4,cycle_count) = Oxidation_Summary(3,cycle_count)/propylene_mass_separated(cycle_count)/(1e6); % MJ/kg

    Reduction_Summary(1,cycle_count) = round(mean(Reduction_Potential_temp),1); % Volt
    Reduction_Summary(2,cycle_count) = abs(trapz(CA_time(Reduction_index_temp(1):Reduction_index_temp(2)),CA_Current_Smooth(Reduction_index_temp(1):Reduction_index_temp(2)))/1000); % Coulomb
    Reduction_Summary(3,cycle_count) = abs(Reduction_Summary(2,cycle_count)*Reduction_Summary(1,cycle_count)); % Joules
    Reduction_Summary(4,cycle_count) = abs(Reduction_Summary(3,cycle_count)/propylene_mass_separated(cycle_count))/(1e6); % MJ/kg

    cycle_count = cycle_count + 1;
end

end