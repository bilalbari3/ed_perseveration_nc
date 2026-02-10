function process_data_table(tbl_batch, save_name)

%% process data at individual provider level
all_providers = unique(tbl_batch.Provider); % enumerate all providers
n_provider = length(all_providers); % how many unique providers
n_unique_order = size(tbl_batch.Order{1}, 2); % cardinality of order set
order_tempate = zeros(1, n_unique_order); % preallocate 

h_shift = 12; % hours to define a shift (NaN first order at start of shift)
h_order_volume = 0.5; % hours to define order/provider count window
m_order_per_min = 5; % how many orders by current provider in past X minutes

% load marginal order distribution
load('/Users/bib002/Dropbox/research/data/ed_perseveration/KL_div_data.mat','Po')
Po_diag = diag(Po); % precompute diagonal of marginal order distribution to generate chance levels of perseveration

% preallocate all analysis variables
persev = NaN(height(tbl_batch), 1);
persev_chance = NaN(height(tbl_batch), 1);
census = NaN(height(tbl_batch), 1);
patient = NaN(height(tbl_batch), 1);
age = NaN(height(tbl_batch), 1);
sex = categorical(NaN(height(tbl_batch), 1));
race = categorical(NaN(height(tbl_batch), 1));
ethnicity = categorical(NaN(height(tbl_batch), 1));
language = categorical(NaN(height(tbl_batch), 1));
provider = NaN(height(tbl_batch), 1);
order_time = NaT(height(tbl_batch), 1);
interorder_interval = NaN(height(tbl_batch), 1);
n_orders = NaN(height(tbl_batch), 1);
n_orders_prev = NaN(height(tbl_batch), 1);
provider_volume = NaN(height(tbl_batch), 1);
order_volume = NaN(height(tbl_batch), 1);
order_per_time = NaN(height(tbl_batch), 1);

% initialize indexing
persev_ind = 1; 
for curr_provider_ind = 1:n_provider
    fprintf('On %i of %i\n', curr_provider_ind, n_provider)
    curr_provider = all_providers(curr_provider_ind);
    tbl_provider = tbl_batch(tbl_batch.Provider == curr_provider, :);

    curr_orders = cell2mat(tbl_provider.Order);
    if size(curr_orders, 1) > 1 % if just 1 order for this provider, skip it and leave persev as NaN'd
        tmp_persev = sum(curr_orders(1:end-1, :) .* curr_orders(2:end, :), 2); % define persev as the inner product
        tmp_persev = [NaN; tmp_persev];

        % assign all easily assigned variables
        indices = persev_ind:persev_ind + height(tbl_provider) - 1;
        persev(indices) = tmp_persev;
        census(indices) = tbl_provider.Census;
        patient(indices) = tbl_provider.Patient;
        age(indices) = tbl_provider.Age;
        sex(indices) = tbl_provider.Sex;
        race(indices) = tbl_provider.Race;
        ethnicity(indices) = tbl_provider.Ethnicity;
        language(indices) = tbl_provider.Language;
        provider(indices) = tbl_provider.Provider;
        order_time(indices) = tbl_provider.Time;
        interorder_interval(indices) = [NaN; minutes(diff(tbl_provider.Time))];
        n_orders(indices) = sum(curr_orders, 2);
        n_orders_prev(indices) = [NaN; sum(curr_orders(1:end - 1, :), 2)];

        % expected perseveration by chance
        P_nooverlap = 1 - sum(curr_orders * Po_diag, 2); % probability of no perseveration
        P_nooverlap = [NaN; P_nooverlap(1:end - 1)]; % bump over by one trial
        P_overlap = 1 - P_nooverlap.^sum(curr_orders, 2); % probability of at least 1 perseveration
        persev_chance(indices) = P_overlap;


        tmp_provider_volume = [];
        tmp_order_volume = [];
        tmp_order_per_time = [];
        for j = 1:height(tbl_provider)
            if mod(j, 100) == 0
                fprintf('\tOn %i of %i\n', j, height(tbl_provider))
            end
            curr_order_time = tbl_provider.Time(j);

            % provider and order volume on the basis of all providers
            provider_count_cutoff_mask = tbl_batch.Time > curr_order_time - hours(h_order_volume) & tbl_batch.Time <= curr_order_time; % use h_order_volume hours to decide on number of providers in ED
            tmp_provider_volume = [tmp_provider_volume; numel(unique(tbl_batch.Provider(provider_count_cutoff_mask)))]; % number of unique providers

            tmp_order_volume = [tmp_order_volume; sum(cell2mat(tbl_batch.Order(provider_count_cutoff_mask)), 'all')];

            % order per time on the basis of current provider
            mask_order_per_time = tbl_provider.Time > curr_order_time - minutes(m_order_per_min) & ...
                                  tbl_provider.Time <= curr_order_time;
            tmp_order_per_time = [tmp_order_per_time; ...
                                  sum(cell2mat(tbl_provider.Order(mask_order_per_time)), 'all')];
        end
        provider_volume(indices) = tmp_provider_volume;
        order_volume(indices) = tmp_order_volume;
        order_per_time(indices) = tmp_order_per_time;

    end
    persev_ind = persev_ind + height(tbl_provider);
end

data_table = table(persev, persev_chance, census, ...
                   categorical(patient), age, sex, race, ethnicity, language, ...
                   categorical(provider), ...
                   order_time, interorder_interval, ...
                   n_orders, n_orders_prev, ...
                   provider_volume, order_volume, order_per_time, ...
                   'VariableNames', ...
                  {'Persev','Persev_Chance','Census',...
                   'Patient','Age','Sex','Race','Ethnicity','Language', ...
                   'Provider',...
                   'Order_Time','Interorder_Interval','N_Orders','N_Orders_Prev', ...
                   'Provider_Volume','Order_Volume','Order_Per_Time'});

save(fullfile('/Users/bib002/Dropbox/research/data/ed_perseveration/', save_name), 'data_table')
fprintf('Saved data_table\n')