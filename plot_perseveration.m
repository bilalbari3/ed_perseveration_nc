%% 
load('/Users/bib002/Dropbox/research/data/ed_perseveration/data_table_all_order_number.mat')
dt = data_table; clear data_table
persev_binary = dt.Persev; persev_binary(persev_binary >= 1) = 1;
dt.("Persev_Binary") = persev_binary;
dt.("Order_Per_Patient") = dt.Order_Volume ./ dt.Census;
dt.("Patient_Per_Provider") = dt.Census ./ dt.Provider_Volume;
dt.("Census_Log") = log(dt.Census);
myColors = importColors_bb();
%% quick analyses
figure; hold on
m_IOI = dt.Interorder_Interval > 16*60;

mask_first_order = double(dt.Patient_First_Order == 1);
% mask_other_order = double(dt.Order_Rel_First_BACKWARDS == -1);
mask_other_order = double(dt.Order_Rel_First_FORWARDS == 1);

y = dt.Persev_Binary;

mask_other_order(m_IOI) = NaN;
mask_first_order(m_IOI) = NaN;
y(m_IOI) = NaN;

y_other = y(mask_other_order == 1);
y_first = y(mask_first_order == 1);

errorbar(1, mean(y_other, "omitnan"), sem(y_other), 'linewidth', 2, 'Color', myColors.orange)
errorbar(2, mean(y_first, "omitnan"), sem(y_first), 'linewidth', 2, 'Color', myColors.blue)

% mu_y_other = mean(y_other, "omitnan");
% [~, ci_y_other] = binofit(sum(y_other == 1), sum(~isnan(y_other)));
% errorbar(1, mu_y_other, mu_y_other - ci_y_other(1), ci_y_other(2) - mu_y_other)
% 
% mu_y_first = mean(y_first, "omitnan");
% [~, ci_y_first] = binofit(sum(y_first == 1), sum(~isnan(y_first)));
% errorbar(2, mu_y_first, mu_y_first - ci_y_first(1), ci_y_first(2) - mu_y_first)

xlim([0 3])
ylim([0.05 0.3])

%%

bin_edges = -2:3;
for bin_index = 1:length(bin_edges) - 1
    bin_left = bin_edges(bin_index);
    bin_right = bin_edges(bin_index + 1);
    if bin_left <= 0
        x = x_back;
    else
        x = x_for;
    end

    mask_x = x >= bin_left & x < bin_right;
    y_plot = y(mask_x);
    errorbar(bin_left, mean(y_plot, "omitnan"), sem(y_plot), ...
        'linewidth', 2, 'Color', myColors.black)
end

xlabel('Order relative to first')
ylabel('Persev')

%%
figure; hold on
m_IOI = dt.Interorder_Interval > 16*60;
m_shift_start = dt.Time_In_Shift == 0;

relative_trials = [-3 2];

ind_back = strfind(dt.Order_Rel_First_BACKWARDS', relative_trials(1):0);
ind_for = strfind(dt.Order_Rel_First_FORWARDS', 0:relative_trials(2));

ind_intersect = intersect(ind_back - relative_trials(1), ind_for);

for trial_index = relative_trials(1):relative_trials(2)
    last_orders = dt.Order_Rel_Last_FORWARDS(ind_intersect + trial_index);
    ind_intersect = ind_intersect(last_orders ~= 0); % no last orders placed
end

y = dt.Persev_Binary;
y(m_IOI | m_shift_start) = NaN;
for trial_index = relative_trials(1):relative_trials(2)
    y_plot = y(ind_intersect + trial_index);
    errorbar(trial_index, mean(y_plot, "omitnan"), sem(y_plot), ...
        'linewidth', 2, 'Color', myColors.black)
end

%%
figure; hold on
m_IOI = dt.Interorder_Interval <= 16*60;
% m_last = dt.Patient_Last_Order == 0;

x_back = dt.Order_Rel_First_BACKWARDS;
x_for = dt.Order_Rel_First_FORWARDS;
y = dt.Persev_Binary;

x_back = x_back(m_IOI);% & m_last);
x_for = x_for(m_IOI);% & m_last);
y = y(m_IOI);% & m_last);

bin_edges = -2:3;
for bin_index = 1:length(bin_edges) - 1
    bin_left = bin_edges(bin_index);
    bin_right = bin_edges(bin_index + 1);
    if bin_left <= 0
        x = x_back;
    else
        x = x_for;
    end

    mask_x = x >= bin_left & x < bin_right;
    y_plot = y(mask_x);
    errorbar(bin_left, mean(y_plot, "omitnan"), sem(y_plot), ...
        'linewidth', 2, 'Color', myColors.black)
end

xlabel('Order relative to first')
ylabel('Persev')