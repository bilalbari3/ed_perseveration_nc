%% load data
load('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_all.mat')
load('/Users/bib002/Dropbox/research/data/ed_perseveration/data_table_all.mat')

%%
all_patients = unique(tbl_batch.Patient);
n_patients = length(all_patients);

n_data_table = height(data_table);
first_order = zeros(n_data_table, 1);
last_order = zeros(n_data_table, 1);

for patient_ind = 1:n_patients
    fprintf('On patient %i of %i\n', patient_ind, n_patients)
    curr_patient = all_patients(patient_ind);
    tbl_patient = tbl_batch(tbl_batch.Patient == curr_patient, :);

    % find first order for patient
    mask_first_order = data_table.Patient == categorical(tbl_patient.Patient(1)) & ... % patient must match
                       data_table.Provider == categorical(tbl_patient.Provider(1)) & ... %% provider must match
                       data_table.Order_Time == tbl_patient.Time(1);

    % find last order for patient
    mask_last_order = data_table.Patient == categorical(tbl_patient.Patient(end)) & ... % patient must match
                       data_table.Provider == categorical(tbl_patient.Provider(end)) & ... %% provider must match
                       data_table.Order_Time == tbl_patient.Time(end);

    if height(tbl_patient) == 1 % if patient only has 1 order, NaN because otherwise first and last order are identical
        first_order(mask_first_order) = NaN;    
        last_order(mask_last_order) = NaN;
    else
        first_order(mask_first_order) = 1;    
        last_order(mask_last_order) = 1;
    end
end

%%
data_table.("Patient_First_Order") = first_order;
data_table.("Patient_Last_Order") = last_order;

%%
save('/Users/bib002/Dropbox/research/data/ed_perseveration/data_table_all_firstlastorder.mat', 'data_table')


%% some rudimentary plots
x_orig = dt.Interorder_Interval;
y_orig = dt.Persev_Binary;

mask_IOI = dt.Interorder_Interval > 16*60;
mask_IOI = ~mask_IOI;
mask_first1 = dt.Patient_First_Order == 1; % first order only
mask_first0 = dt.Patient_First_Order == 0; % non-first orders
diff_patient = [NaN; diff(double(dt.Patient))];
mask_persev_same_patient = diff_patient == 0 & dt.Persev_Binary == 1; % if perseveration was to _same_ patient
mask_persev_diff_patient = diff_patient == 1 & dt.Persev_Binary == 1; % if perseveration was to _different_ patient

figure; hold on
bin_edges = 1:1:60; 

% x = x_orig; y = y_orig;
% x(mask_persev_same_patient) = NaN; y_mask_persev_same_patient = NaN;
% x = x(mask_IOI);
% y = y(mask_IOI);
x = x_orig; x = x(mask_IOI & mask_first1);
y = y_orig; y = y(mask_IOI & mask_first1);
for bin_index = 1:length(bin_edges) - 1
    bin_left = bin_edges(bin_index);
    bin_right = bin_edges(bin_index + 1);
    y_ebar = y(x >= bin_left & x < bin_right);
    errorbar(mean([bin_left bin_right], "omitnan"), mean(y_ebar, "omitnan"), sem(y_ebar), ...
        'linewidth', 2, 'Color', myColors.blue)
end

% x = x_orig; y = y_orig;
% x(mask_persev_diff_patient) = NaN; y_mask_persev_diff_patient = NaN;
% x = x(mask_IOI);
% y = y(mask_IOI);
x = x_orig; x = x(mask_IOI & mask_first0);
y = y_orig; y = y(mask_IOI & mask_first0);
for bin_index = 1:length(bin_edges) - 1
    bin_left = bin_edges(bin_index);
    bin_right = bin_edges(bin_index + 1);
    y_ebar = y(x >= bin_left & x < bin_right);
    errorbar(mean([bin_left bin_right], "omitnan"), mean(y_ebar, "omitnan"), sem(y_ebar), ...
        'linewidth', 2, 'Color', myColors.orange)
end