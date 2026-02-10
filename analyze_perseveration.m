clear; close all; clc
load('/Users/bib002/Dropbox/research/data/ed_perseveration/data_table_all_order_number.mat')
dt = data_table; clear data_table
fprintf('Loaded data table\n')

myColors = importColors_bb();
saveFigLoc = '/Users/bib002/Documents/git-repositories/ed_perseveration/figures';
saveDataLoc = '/Users/bib002/Dropbox/research/data/ed_perseveration/';

%% preanalyze some variables
persev_binary = dt.Persev; persev_binary(persev_binary >= 1) = 1;
dt.("Persev_Binary") = persev_binary;
dt.("Order_Per_Patient") = dt.Order_Volume ./ dt.Census;
dt.("Patient_Per_Provider") = dt.Census ./ dt.Provider_Volume;
dt.("Census_Log") = log(dt.Census);


%% plot it
% default: NaN out ioi > 12 hours
% persev: binary, fraction
% covariates: number of orders in batch, baseline probability
% control: last and current trial batch size = 1
% control: NaN out perseveration to same patient

% cog variables: Census, Order_Volume, Order_Per_Patient, Patient_Per_Provider, Provider_Volume
% time variables: Interorder_Interval, Order_Per_Time

% cog_load_variable = 'Census';
% cog_load_variable = 'Census_Log';
% cog_load_variable = 'Interorder_Interval';
% cog_load_variable = 'Order_Per_Time';
% cog_load_variable = 'Time_In_Shift';
cog_load_variable = 'Patient_First_Order';

fit_model_flag = true;
[f_handle, mod, sp_handle] = plot_ed_persev(dt, cog_load_variable, ...
                                 'Color', myColors.bluishGreen, ...
                                 'Fit_Model_Flag', fit_model_flag, ...
                                 'Model_Class', 'mixed');
% saveFigureIteration(f_handle, saveFigLoc, cog_load_variable, 'AutoSave_Flag', false)
if fit_model_flag == true
    save(fullfile(saveDataLoc, 'regressions', strcat(cog_load_variable, '_regressions.mat')), 'mod','-v7.3')
end

%% run stats
mask_ioi = dt.Interorder_Interval > 16*60;
mask_first_trial = dt.Time_In_Shift == 0;
diff_patient = [NaN; diff(double(dt.Patient))];
mask_persev_same_patient = diff_patient == 0 & dt.Persev_Binary == 1;
mask_orders = dt.N_Orders_Prev == 1 & dt.N_Orders == 1; % only if prior order set and current order set = 1;
mask_orders = ~mask_orders; % flip mask to NaN out instances where order set size is not 1

mask_global = mask_ioi | mask_first_trial;
% mask_global = mask_ioi | mask_persev_same_patient;
% mask_global = mask_ioi | mask_persev_same_patient | mask_orders;

y = dt.Persev_Binary;
y(mask_global) = NaN;

dt_short = table(y, dt.Census, dt.Interorder_Interval, dt.Time_In_Shift, ...
                 dt.Patient_First_Order, dt.Persev_Chance, ...
                 dt.Provider, ...
                 'VariableNames', ...
                {'Persev','Census','IOI','Time','First_Order', 'Persev_Chance', 'Provider'});
dt_short.Census(mask_global) = NaN;
dt_short.IOI(mask_global) = NaN;
dt_short.Time(mask_global) = NaN;
dt_short.First_Order(mask_global) = NaN;
dt_short.Persev_Chance(mask_global) = NaN;

dt_short.Census = normalize(dt_short.Census);
dt_short.IOI = normalize(dt_short.IOI);
dt_short.Time = normalize(dt_short.Time);
dt_short.First_Order = normalize(dt_short.First_Order);
dt_short.Persev_Chance = normalize(dt_short.Persev_Chance);

%%
modglobal = fitglme(dt_short, 'Persev ~ Census + IOI + Time + First_Order + Persev_Chance + (1|Provider)', ...
                    'distribution','binomial');
save(fullfile(saveDataLoc, 'regressions', 'Global_regressions.mat'), 'modglobal','-v7.3')

%%
function [f_handle, mod, t] = plot_ed_persev(dt, x_var, varargin)
    p = inputParser;
    p.addParameter('Fit_Model_Flag', false)
    p.addParameter('Model_Class', 'fixed') % fixed or mixed
    p.addParameter('Color', [0 0 0]);
    p.parse(varargin{:});

    f_handle = figure('Position',[90 303 1478 653]);
    sgtitle(x_var,'fontsize',18,'interpreter','none')
    n_bins = 11;
    i_remap = [1 4 2 5 3 6]; % remap subplots to make more sense
    x_plot = dt.(x_var);

    % interorder-interval mask
    mask_ioi = dt.Interorder_Interval > 16*60; % more than X hours defines next shift; nan out 1st order
    x_plot(mask_ioi) = NaN;

    % first trial of shift
    mask_first_trial = dt.Time_In_Shift == 0;
    x_plot(mask_first_trial) = NaN;

    % perseverate to same patient mask
    diff_patient = [NaN; diff(double(dt.Patient))];
    mask_persev_same_patient = diff_patient == 0 & dt.Persev_Binary == 1; % if perseveration was to _same_ patient
    % mask_persev_same_patient = diff_patient == 0;

    % order set size 1
    mask_orders = dt.N_Orders_Prev == 1 & dt.N_Orders == 1; % only if prior order set and current order set = 1;
    mask_orders = ~mask_orders; % flip mask to NaN out instances where order set size is not 1
   
    mod = struct();

    % determine plot type
    switch x_var
        case {'Census','Census_Log'}
            plot_type = 'binned';
        case {'Interorder_Interval','Order_Per_Time','Time_In_Shift'}
            plot_type = 'errorbar';
        case {'Patient_First_Order'}
            plot_type = 'binary';
    end

    for i = 1:6
        switch i
            case 1
                title_label = 'Persev_Binary';
                y_plot = dt.Persev_Binary;
                y_plot(mask_ioi) = NaN;
                model_type = 'logistic';
            case 2
                title_label = 'Persev_Fraction';
                y_plot = dt.Persev./dt.N_Orders_Prev;
                y_plot(mask_ioi) = NaN;
                model_type = 'linear';
            case 3
                title_label = 'Persev_Binary: Different Patients';
                y_plot = dt.Persev_Binary;
                y_plot(mask_ioi | mask_persev_same_patient) = NaN;
                model_type = 'logistic';
            case 4
                title_label = 'Persev_Fraction: Different Patients';
                y_plot = dt.Persev./dt.N_Orders_Prev;
                y_plot(mask_ioi | mask_persev_same_patient) = NaN;
                model_type = 'linear';
            case 5
                title_label = 'Persev_Binary: Order Set Size 1';
                y_plot = dt.Persev_Binary;
                y_plot(mask_ioi | mask_orders) = NaN;
                model_type = 'logistic';
            case 6
                title_label = sprintf('Persev_Binary: Different Patients\nand Order Set Size 1');
                y_plot = dt.Persev_Binary;
                y_plot(mask_ioi | mask_persev_same_patient | mask_orders) = NaN;
                model_type = 'logistic'; % reduces to logistic under these conditions
        end
        t(i) = subplot(2,3,i_remap(i)); hold on
        title(title_label,'interpreter','none')
        if strcmp(plot_type, 'binned')
            plot_my_binned(y_plot)
        elseif strcmp(plot_type, 'errorbar')
            plot_my_ebar(y_plot)
        elseif strcmp(plot_type, 'binary')
            plot_binary(y_plot)
        end
        if p.Results.Fit_Model_Flag == true
            fprintf('\tFitting model %i\n', i)
            mod_name = strcat('mod', num2str(i));
            mod.(mod_name) = fit_my_model(y_plot, model_type);
        end
    end

    for cp = t % standardize each subplot
        subplot(cp)
        set(cp,'fontsize',18)
        xlabel(x_var,'interpreter','none')
        ylabel('P(repeat order)')
    end

    function plot_my_binned(y_plot)
        % plot percentile binned data
        plotBinned(x_plot, y_plot, n_bins, gca, 'Color', p.Results.Color);
    end

    function plot_my_ebar(y_plot)
        % plot errorbar
        switch x_var
            case {'Interorder_Interval'}
                bin_edges = 1:1:45; 
            case {'Order_Per_Time'}
                bin_edges = 1:30;
            case ('Time_In_Shift')
                bin_edges = 0:0.5:8;
            otherwise
                error('This input not supported')
        end
        for bin_index = 1:length(bin_edges) - 1
            bin_left = bin_edges(bin_index);
            bin_right = bin_edges(bin_index + 1);
            y_ebar = y_plot(x_plot >= bin_left & x_plot < bin_right);
            errorbar(mean([bin_left bin_right], "omitnan"), mean(y_ebar, "omitnan"), sem(y_ebar), ...
                'linewidth', 2, 'Color', p.Results.Color)
        end
    end

    function plot_binary(y_plot) 
        % plot first order for new patient
        mask_first_order = double(dt.Patient_First_Order == 1); % first order for patient
        mask_other_order = double(dt.Patient_First_Order == 0); % immediate preceding order
        y_other = y_plot(mask_other_order == 1);
        y_first = y_plot(mask_first_order == 1);

        mu_y_other = mean(y_other, "omitnan");
        % [~, ci_y_other] = binofit(sum(y_other == 1), sum(~isnan(y_other)));
        bar(0, mu_y_other, 'FaceColor', p.Results.Color, 'EdgeColor', 'none')
        % errorbar(0, mu_y_other, mu_y_other - ci_y_other(1), ci_y_other(2) - mu_y_other, ...
        %     'linewidth', 2, 'Color', [0 0 0])
        errorbar(0, mu_y_other, sem(y_other), ...
            'linewidth', 2, 'Color', [0 0 0])

        mu_y_first = mean(y_first, "omitnan");
        % [~, ci_y_first] = binofit(sum(y_first == 1), sum(~isnan(y_first)));
        bar(1, mu_y_first, 'FaceColor', p.Results.Color, 'EdgeColor', 'none')
        % errorbar(1, mu_y_first, mu_y_first - ci_y_first(1), ci_y_first(2) - mu_y_first, ...
        %     'linewidth', 2, 'Color', [0 0 0])
        errorbar(1, mu_y_first, sem(y_first), ...
            'linewidth', 2, 'Color', [0 0 0])

        xlim([-1 2])
    end
    
    function mod_out = fit_my_model(y_plot, model_type)
        switch model_type
            case 'logistic'
                % do not normalize y_plot
                dt_analysis = table(normalize(x_plot), y_plot, normalize(dt.Persev_Chance), dt.Provider, ... 
                                    'VariableNames',...
                                    {x_var, 'Persev', 'Persev_Chance', 'Provider'});

                if strcmp(p.Results.Model_Class, 'fixed')
                    mod_out = fitglm(dt_analysis, sprintf('Persev ~ %s + Persev_Chance', x_var), ...
                                    'Distribution', 'binomial');
                elseif strcmp(p.Results.Model_Class, 'mixed')
                    mod_out = fitglme(dt_analysis, sprintf('Persev ~ %s + Persev_Chance + (1|Provider)', x_var), ...
                                  'Distribution', 'binomial');
                end
            case 'linear'
                % normalize y_plot
                dt_analysis = table(normalize(x_plot), normalize(y_plot), normalize(dt.Persev_Chance), dt.Provider, ... 
                                    'VariableNames',...
                                    {x_var, 'Persev', 'Persev_Chance', 'Provider'});
                if strcmp(p.Results.Model_Class, 'fixed')
                    mod_out = fitlm(dt_analysis, sprintf('Persev ~ %s + Persev_Chance', x_var));
                elseif strcmp(p.Results.Model_Class, 'mixed')
                    mod_out = fitlme(dt_analysis, sprintf('Persev ~ %s + Persev_Chance + (1|Provider)', x_var));
                end
        end
    end
end

%% regressions vs time
time = dt.Order_Time;
time = dateshift(time, 'start', 'week'); % strip hours and minutes

census = dt.Census;
persev = dt.Persev; persev(persev >= 1) = 1;
persev_chance = dt.Persev_Chance;
provider = dt.Provider;

[time_sorted, sort_inds] = sort(unique(time));
time_sorted(isnat(time_sorted)) = []; % drop NaT values
coef = [];
for curr_date = time_sorted'
    fprintf('On %s\n', curr_date)
    mask_date = time == curr_date;
    c_census = census(mask_date);
    c_persev = persev(mask_date);
    c_persev_chance = persev_chance(mask_date);
    c_provider = provider(mask_date);
    dt_short = table(c_census, c_persev, c_persev_chance, c_provider, ...
        'VariableNames',{'Census','Persev','Persev_Chance','Provider'});
    % mod = fitglm(c_census, c_persev, 'distribution', 'binomial');
    mod = fitglm(dt_short, 'Persev ~ Census','distribution','binomial');
    % mod = fitglme(dt_short, 'Persev ~ Census + Persev_Chance + (1|Provider)','distribution','binomial');
    coef = [coef; mod.Coefficients.Estimate(2)];
end

plot(time_sorted, coef, 'x')

%% analyze using all interesting variables
Persev_Binary = dt.Persev_Binary;
Census = dt.Census;
Census_Log = dt.Census_Log;
Interorder_Interval = dt.Interorder_Interval;
Order_Per_Time = dt.Order_Per_Time;
Persev_Chance = dt.Persev_Chance;
Provider = dt.Provider;

% mask based on these variables
Max_Time_In_Shift = dt.Max_Time_In_Shift;
Time_In_Shift = dt.Time_In_Shift;

mask_nan = Time_In_Shift == 0 | Interorder_Interval > 16*60; % get rid of first order of shift, and first order placed after 16 hours
Persev_Binary(mask_nan) = NaN;

dt_mod = table(Persev_Binary, normalize(Census), normalize(Census_Log), normalize(Interorder_Interval), normalize(Order_Per_Time), ...
                    normalize(Time_In_Shift), normalize(Persev_Chance), categorical(Provider), ...
                    'VariableNames', ...
                    {'Persev_Binary','Census', 'Census_Log','Interorder_Interval','Order_Per_Time', ...
                     'Time_In_Shift','Persev_Chance','Provider'});

mod_global = fitglm(dt_mod, ['Persev_Binary ~ Census + Interorder_Interval + Order_Per_Time + Time_In_Shift + ' ...
    'Persev_Chance'], ...
       'distribution','binomial')