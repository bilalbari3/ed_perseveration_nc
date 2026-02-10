clear; close all; clc
load('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_all.mat')
fprintf('Loaded data\n')

%% generate order matrix for each patient — calculate first N orders
% this is to calculate KL-divergence and H(order | patient)
rng(1) % because we randomly remove orders above threshold
all_patients = unique(tbl_batch.Patient)';
n_patients = length(all_patients);

n_unique_order = size(tbl_batch.Order{1}, 2);

order_matrix = NaN(n_patients, n_unique_order); % initialize order matrix as NaN
census = NaN(n_patients, 1); % average census at time of each order
time = NaT(n_patients, 1); % average datetime at time of each order
patient = categorical(NaN(n_patients, 1)); % patient ID
provider = cell(n_patients, 1);
age = NaN(n_patients, 1);
sex = categorical(NaN(n_patients, 1));
race = categorical(NaN(n_patients, 1));
ethnicity = categorical(NaN(n_patients, 1));
language = categorical(NaN(n_patients, 1));
number_ordering_providers = NaN(n_patients, 1); % number of ordering providers for first N_min_orders orders

N_min_orders = 10;

for curr_patient_ind = 1:n_patients
    fprintf('On %i of %i\n', curr_patient_ind, n_patients)
    curr_patient = all_patients(curr_patient_ind);
    tbl_patient = tbl_batch(tbl_batch.Patient == curr_patient, :); % table for just that patient
    curr_order_matrix = cell2mat(tbl_patient.Order);

    curr_total_orders = sum(curr_order_matrix, 'all'); % how many orders for this patient?
    if curr_total_orders >= N_min_orders % at least N_min_orders for this patient
        index_order_threshold = 1; % set this to 1 for census calculation
        if curr_total_orders > N_min_orders % if more than N_min_orders, we need to subsample
            curr_order_cumsum = cumsum(sum(curr_order_matrix, 2)); % cumulative sum
            index_order_threshold = find(curr_order_cumsum  >= N_min_orders, 1, 'first');
            order_threshold = find(curr_order_matrix(index_order_threshold, :)); % what orders were made for threshold set
            n_excess = curr_order_cumsum(index_order_threshold) - N_min_orders; % how many orders to discard?
            
            order_threshold_subsamp = randsample(order_threshold, length(order_threshold) - n_excess); % remove excess orders
            order_vector = zeros(1, n_unique_order);
            order_vector(order_threshold_subsamp) = 1; % order vector for that last row
            curr_order_matrix(index_order_threshold, :) = order_vector;
            curr_order_matrix = curr_order_matrix(1:index_order_threshold, :);
        end
        order_matrix(curr_patient_ind, :) = sum(curr_order_matrix, 1);
        census(curr_patient_ind) = mean(tbl_patient.Census(1:index_order_threshold));
        time(curr_patient_ind) = mean(tbl_patient.Time(1:index_order_threshold));
        patient(curr_patient_ind) = categorical(tbl_patient.Patient(1));
        age(curr_patient_ind) = tbl_patient.Age(1);
        sex(curr_patient_ind) = tbl_patient.Sex(1);
        race(curr_patient_ind) = tbl_patient.Race(1);
        ethnicity(curr_patient_ind) = tbl_patient.Ethnicity(1);
        language(curr_patient_ind) = tbl_patient.Language(1);

        ordering_providers = tbl_patient.Provider(1:index_order_threshold); % pull out providers that placed orders
        ordering_providers = unique(ordering_providers); % who were the providers?
        number_ordering_providers(curr_patient_ind) = length(ordering_providers); % append number of providers
        provider{curr_patient_ind} = ordering_providers;
    end
end

% generate order matrix for order categories
n_unique_order_cat = size(tbl_batch.Order_Category{1}, 2);
order_cat_matrix = NaN(n_patients, n_unique_order_cat); % initialize order category matrix as NaN
load('/Users/bib002/Dropbox/research/data/ed_perseveration/order_to_cat_LUT.mat') % lookup table for order to category

for curr_patient_ind = 1:n_patients % go through each patient
    fprintf('On %i of %i\n', curr_patient_ind, n_patients)
    if ~isnan(census(curr_patient_ind)) % if this is not a nan entry
        order_inds = find(order_matrix(curr_patient_ind, :)); % which unique orders were placed
        order_cat_inds_unique = order_to_cat_LUT(order_inds); % convert orders to order categories
        [order_cat_sum, order_cat_inds] = groupcounts(order_cat_inds_unique'); % how many entries for each category?
        order_cat_vector = zeros(1, n_unique_order_cat);
        order_cat_vector(order_cat_inds) = order_cat_sum;
        order_cat_matrix(curr_patient_ind, :) = order_cat_vector;
    end
end

% order matrix for order categories 2
order_cat2_matrix = NaN(n_patients, 2);
order_cat2_matrix(:, 1) = order_cat_matrix(:, 33);
order_cat2_matrix(:, 2) = sum(order_cat_matrix(:, [1:32 34:end]), 2);

%% analyze KL divergence and H(orders | patient) for individual orders, order categories, and order categories 2
% all orders
cond_ent = NaN(length(census), 1);
KL_div = NaN(length(census), 1);
MI = NaN(length(census), 1);
Po = mean(order_matrix, "omitnan"); Po = Po./sum(Po);
H_Po = -1*sum(Po.*log10(Po),'omitnan');
for n = 1:length(cond_ent)
    if ~isnan(census(n)) % if not a nan entry
        tmp_orders = order_matrix(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders); % convert to probabilities
        cond_ent(n) = -1*sum(tmp_orders.*log10(tmp_orders), "omitnan");
        KL_div(n) = kldiv(1:n_unique_order, ...
                          tmp_orders + eps, ...
                          Po + eps);
        MI(n) = H_Po - cond_ent(n);
    end
end

% order categories
cond_ent_cat = NaN(length(census), 1);
KL_div_cat = NaN(length(census), 1);
MI_cat = NaN(length(census), 1);
Po_cat = mean(order_cat_matrix, "omitnan"); Po_cat = Po_cat./sum(Po_cat);
H_Po_cat = -1*sum(Po_cat.*log10(Po_cat),'omitnan');
for n = 1:length(cond_ent_cat)
    if ~isnan(census(n)) % if not a nan entry
        tmp_orders = order_cat_matrix(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders); % convert to probabilities
        cond_ent_cat(n) = -1*sum(tmp_orders.*log10(tmp_orders), "omitnan");
        KL_div_cat(n) = kldiv(1:n_unique_order_cat, ...
                              tmp_orders + eps, ...
                              Po_cat + eps);
        MI_cat(n) = H_Po_cat - cond_ent_cat(n);
    end
end

% order categories 2
cond_ent_cat2 = NaN(length(census), 1);
KL_div_cat2 = NaN(length(census), 1);
MI_cat2 = NaN(length(census), 1);
Po_cat2 = mean(order_cat2_matrix, "omitnan"); Po_cat2 = Po_cat2./sum(Po_cat2);
H_Po_cat2 = -1*sum(Po_cat2.*log10(Po_cat2),'omitnan');
for n = 1:length(cond_ent_cat2)
    if ~isnan(census(n)) % if not a nan entry
        tmp_orders = order_cat2_matrix(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders); % convert to probabilities
        cond_ent_cat2(n) = -1*sum(tmp_orders.*log10(tmp_orders), "omitnan");
        KL_div_cat2(n) = kldiv(1:2, ...
                               tmp_orders + eps, ...
                               Po_cat2 + eps);
        MI_cat2(n) = H_Po_cat2 - cond_ent_cat2(n);
    end
end

%% save data
save('/Users/bib002/Dropbox/research/data/ed_perseveration/KL_div_data.mat', ...
    '-regexp', '^(?!(tbl_batch)$).','-v7.3') % save except tbl_batch

%% generate data table
dt = table(census, order_matrix, order_cat_matrix, order_cat2_matrix, ...
           cond_ent, KL_div, MI, cond_ent_cat, KL_div_cat, MI_cat, cond_ent_cat2, KL_div_cat2, MI_cat2, ...
           time, patient, age, sex, race, ethnicity, language, number_ordering_providers, provider, ...
    'VariableNames', ...
    {'Census','Order','Order_Category','Order_Category2','Cond_Ent','KL_Div','MI',...
     'Cond_Ent_Cat','KL_Div_Cat','MI_Cat','Cond_Ent_Cat2','KL_Div_Cat2','MI_Cat2',...
     'Time','Patient','Age','Sex','Race','Ethnicity','Language','Number_Ordering_Providers','Providers'});
dt_struct = struct();
dt_struct.Po = Po;
dt_struct.Po_cat = Po_cat;
dt_struct.Po_cat2 = Po_cat2;
dt_struct.n_unique_order = n_unique_order;
dt_struct.n_unique_order_cat = n_unique_order_cat;

save('/Users/bib002/Dropbox/research/data/ed_perseveration/KL_div_data_table.mat','dt','dt_struct','-v7.3')


