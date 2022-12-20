function Graph_Results_Summary_v01

Varying_Ni2D_Loading = importdata("Varying_Ni2D_Loading.mat");
Varying_IL_Loading = importdata("Varying_IL_Loading.mat");
Varying_Feed_Composition = importdata("Varying_Feed_Composition.mat");
Varying_Reductive_Potential = importdata("Varying_Reductive_Potential.mat");
Varying_Oxidative_Potential = importdata("Varying_Oxidative_Potential.mat");
Varying_Feed_Exposure_Time = importdata("Varying_Feed_Exposure_Time.mat");

primary_path_address = "C:\Users\Halcy\Dropbox\Nickel_Dithiolene_Research_Experimental_Data\Data_Trends_Analysis\Deconvoluted_ECM_plots\Deconvoluted_ECM_";
scaler = 1e6; % convert mol to µmol
x_scaler = 0.2;
y_cap = 0.2;

x_range = [-2000 2000];
y_range = [0 40];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Varying Ni2D Loading
error_bar_grapher(Varying_Ni2D_Loading,'[N(Et)]_2[Ni(mnt)_2] Loading [mg/cm^2]','Propylene Separated [µmol/cm^2]',scaler,x_scaler,y_cap);

Data_ID_numbers = [1 2 3 4 5];
Data_Label = [' 0 mg/cm^2';' 2 mg/cm^2';' 3 mg/cm^2';' 4 mg/cm^2';' 8 mg/cm^2'];
Data_Legend_Title = '[N(Et)]_2[Ni(mnt)_2] Loading';

Combined_Grapher(Data_ID_numbers,primary_path_address,'Time Relative to Polarity Switch [s]','Concentration of Separated Propylene [ppm]',...
    x_range,y_range,Data_Label,Data_Legend_Title,0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Varying IL Loading
error_bar_grapher(Varying_IL_Loading,'IL Loading [mg/cm^2]','Propylene Separated [µmol/cm^2]',scaler,x_scaler,y_cap);

Data_ID_numbers = [6 7 5 9 10];
Data_Label = [' 0 mg/cm^2';' 2 mg/cm^2';' 4 mg/cm^2';' 8 mg/cm^2';'16 mg/cm^2'];
Data_Legend_Title = 'Ionic Liquid Loading';

Combined_Grapher(Data_ID_numbers,primary_path_address,'Time Relative to Polarity Switch [s]','Concentration of Separated Propylene [ppm]',...
    x_range,y_range,Data_Label,Data_Legend_Title,4,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Varying Feed Composition
error_bar_grapher(Varying_Feed_Composition,'Propylene Fraction of Feed','Propylene Separated [µmol/cm^2]',scaler,x_scaler,y_cap);

Data_ID_numbers = [1 12 5 14 15];
Data_Label = [' 0.00';' 0.25';' 0.50';' 0.75';' 1.00'];
Data_Legend_Title = 'Propylene Fraction';

Combined_Grapher(Data_ID_numbers,primary_path_address,'Time Relative to Polarity Switch [s]','Concentration of Separated Propylene [ppm]',...
    x_range,y_range,Data_Label,Data_Legend_Title,0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Varying Reductive Potential
error_bar_grapher(Varying_Reductive_Potential,'Reductive Potential vs. Counterside [V]','Propylene Separated [µmol/cm^2]',scaler,1,y_cap);
hold on;
plot([-1.46 -1.46],[0 300],'-.g','linewidth',3,'markersize',30);
xlim([-3.5 0.5]);

Data_ID_numbers = [5 20 1 35];
Data_Label = ['-2.0V';'-3.0V';'-2.0V';'-3.0V'];
Data_Legend_Title = ({'Reductive Potential';
    'Propylene    Propane'});

Combined_Grapher2(Data_ID_numbers,primary_path_address,'Time Relative to Polarity Switch [s]','Gas Concentration [ppm]',...
    x_range,[0 500],Data_Label,Data_Legend_Title,3,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Varying Oxidative Potential
error_bar_grapher(Varying_Oxidative_Potential,'Oxidative Potential vs. Counterside [V]','Propylene Separated [µmol/cm^2]',scaler,x_scaler,y_cap);
hold on;
plot([0.82 0.82],[0 300],'-.m','linewidth',3,'markersize',30);

Data_ID_numbers = [1 1 23 24 25 26 27];
Data_Label = [' 0.0V';'+0.5V';'+1.0V';'+1.5V';'+2.0V';'+2.5V';'+3.0V'];
Data_Legend_Title = ('Oxidative Potential');

Combined_Grapher(Data_ID_numbers,primary_path_address,'Time Relative to Polarity Switch [s]','Concentration of Separated Propylene [ppm]',...
    x_range,y_range,Data_Label,Data_Legend_Title,5,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Varying Feed Exposure

regular_grapher(Varying_Feed_Exposure_Time,'Feed Exposure Time [min]','Propylene Separated [µmol/cm^2]',scaler,x_scaler,y_cap);

Data_ID_numbers = [28 29 30 31 25 33 34];
Data_Label = [' 5 minutes';'10 minutes';'15 minutes';'20 minutes';'30 minutes';'45 minutes';'60 minutes'];
Data_Legend_Title = ('Feed Exposure Time');

Combined_Grapher(Data_ID_numbers,primary_path_address,'Time Relative to Polarity Switch [s]','Concentration of Separated Propylene [ppm]',...
    x_range,y_range,Data_Label,Data_Legend_Title,0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function Combined_Grapher2(ID_vector,primary_path_address,xlabel_name,ylabel_name,x_range,y_range,data_label,legend_title,color_change_index,number_col)

number_plot_lines = length(ID_vector);
gradient_number1 = linspace(1,0,number_plot_lines);

color_blue = [0 0 1];
color_red = [1 0 0];

figure;
hold on;
count = 1;
while count <= number_plot_lines
    if ID_vector(count) < 10
        path_address = strcat(primary_path_address,"0",num2str(ID_vector(count)),".xls");
    else
        path_address = strcat(primary_path_address,num2str(ID_vector(count)),".xls");
    end
    
    color_blue(2) = gradient_number1(count);
    color_red(2) = gradient_number1(count);
    
    graph_plot_array = importdata(path_address);
    if color_change_index == 0
        plot(graph_plot_array(:,1)+120,graph_plot_array(:,2),'color',color_blue,'linewidth',3,'markersize',30);
    else
        if count >= color_change_index
            plot(graph_plot_array(:,1)+120,graph_plot_array(:,2),'color',color_red,'linewidth',3,'markersize',30);
        else
            plot(graph_plot_array(:,1)+120,graph_plot_array(:,2),'color',color_blue,'linewidth',3,'markersize',30);
        end
    end
    count = count + 1;
end
xlim(x_range);
ylim(y_range);
legend(data_label,'NumColumns',number_col,'location','northeast','FontSize',20);
leg = legend('show');
title(leg, legend_title);
ylabel(ylabel_name);
xlabel(xlabel_name);
box on;
set(gca,'FontName','Arial','FontSize',30,'linewidth',2,'TickLength',[0.025 0.025]);
pbaspect([1 1 1]);

end

function Combined_Grapher(ID_vector,primary_path_address,xlabel_name,ylabel_name,x_range,y_range,data_label,legend_title,color_change_index,number_col)

number_plot_lines = length(ID_vector);
gradient_number1 = linspace(1,0,number_plot_lines);

color_blue = [0 0 1];
color_red = [1 0 0];

figure;
hold on;
count = 1;
while count <= number_plot_lines
    if ID_vector(count) < 10
        path_address = strcat(primary_path_address,"0",num2str(ID_vector(count)),".xls");
    else
        path_address = strcat(primary_path_address,num2str(ID_vector(count)),".xls");
    end
    
    color_blue(2) = gradient_number1(count);
    color_red(2) = gradient_number1(count);
    
    graph_plot_array = importdata(path_address);
    if color_change_index == 0
        plot(graph_plot_array(:,1)+120,graph_plot_array(:,2),'color',color_blue,'linewidth',3,'markersize',30);
    else
        if count >= color_change_index
            plot(graph_plot_array(:,1)+120,graph_plot_array(:,2),'color',color_red,'linewidth',3,'markersize',30);
        else
            plot(graph_plot_array(:,1)+120,graph_plot_array(:,2),'color',color_blue,'linewidth',3,'markersize',30);
        end
    end
    count = count + 1;
end
xlim(x_range);
ylim(y_range);
legend(data_label,'NumColumns',number_col,'location','northwest','FontSize',20);
leg = legend('show');
title(leg, legend_title);
ylabel(ylabel_name);
xlabel(xlabel_name);
box on;
set(gca,'FontName','Arial','FontSize',30,'linewidth',2,'TickLength',[0.025 0.025]);
pbaspect([1 1 1]);

end

function regular_grapher(data,xlabel_name,ylabel_name,scaler,xlim_scale,y_range_cap)

number_points = size(data,2);
vector_average = zeros(1,number_points);
vector_standard_error = zeros(1,number_points);
count = 1;
while count <= number_points
    vector_average(count) = mean(data(2:end,count));
    vector_standard_error(count) = std(data(2:end,count))/sqrt(size(data,1)-1);
    count = count + 1;
end

figure;
plot(data(1,:),vector_average*scaler,...
        '.','MarkerSize',40,'MarkerEdgeColor','red','linewidth',5);
xlim([min(data(1,:))-xlim_scale*max(data(1,:)) (1+xlim_scale)*max(data(1,:))]);
ylim([0 y_range_cap]);
ylabel(ylabel_name);
xlabel(xlabel_name);
box on;
set(gca,'FontName','Arial','FontSize',30,'linewidth',2,'TickLength',[0.025 0.025]);
pbaspect([1 1 1]);

end

function error_bar_grapher(data,xlabel_name,ylabel_name,scaler,xlim_scale,y_range_cap)

number_points = size(data,2);
vector_average = zeros(1,number_points);
vector_standard_error = zeros(1,number_points);
count = 1;
while count <= number_points
    vector_average(count) = mean(data(2:end,count));
    vector_standard_error(count) = std(data(2:end,count))/sqrt(size(data,1)-1);
    count = count + 1;
end

figure;
errorbar(data(1,:),vector_average*scaler,vector_standard_error*scaler,...
        '.','MarkerSize',40,'MarkerEdgeColor','red','linewidth',5);
xlim([min(data(1,:))-xlim_scale*max(data(1,:)) (1+xlim_scale)*max(data(1,:))]);
ylim([0 y_range_cap]);
ylabel(ylabel_name);
xlabel(xlabel_name);
box on;
set(gca,'FontName','Arial','FontSize',30,'linewidth',2,'TickLength',[0.025 0.025]);
pbaspect([1 1 1]);

end