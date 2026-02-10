%% process data at individual provider level
all_providers = unique(data_table.Provider); % enumerate all providers
n_provider = length(all_providers); % how many unique providers

n_orders = height(data_table);
order_relative_to_first_backwards = NaN(n_orders, 1);
order_relative_to_first_forwards = NaN(n_orders, 1);
order_relative_to_last_backwards = NaN(n_orders, 1);
order_relative_to_last_forwards = NaN(n_orders, 1);
order_in_shift = NaN(n_orders, 1);
time_in_shift = NaN(n_orders, 1);
max_time_in_shift = NaN(n_orders, 1);
for curr_provider_ind = 1:n_provider
    fprintf('On %i of %i\n', curr_provider_ind, n_provider)
    curr_provider = all_providers(curr_provider_ind);
    curr_provider_indices = find(data_table.Provider == curr_provider);
    tbl_provider = data_table(curr_provider_indices, :);
    n_orders_by_provider = height(tbl_provider);

    if height(tbl_provider) >= 1
    
        tmp_order_relative_to_first_backwards = NaN(n_orders_by_provider, 1);
        tmp_order_relative_to_first_forwards = NaN(n_orders_by_provider, 1);
        tmp_order_relative_to_last_backwards = NaN(n_orders_by_provider, 1);
        tmp_order_relative_to_last_forwards = NaN(n_orders_by_provider, 1);
        tmp_order_in_shift = NaN(n_orders_by_provider, 1);
        tmp_time_in_shift = NaN(n_orders_by_provider, 1);
        tmp_max_time_in_shift = NaN(n_orders_by_provider, 1);
    
        % mask_shift_cutoff = [NaT - NaT; diff(tbl_provider.Order_Time)] >= hours(16); % if IOI > 16 hours, define new shift
        % shift_cutoff_inds = [1; find(mask_shift_cutoff)];
    
        % when the cumulative time in a shift > 16 hours, reset shift
        times = hours(tbl_provider.Order_Time - tbl_provider.Order_Time(1));
        times_of_shift = zeros(size(tbl_provider.Order_Time)); % Initialize the output vector
        threshold = 16;
        cumSum = 0; % Initialize the cumulative sum
        for i = 1:length(times)
            cumSum = cumSum + times(i); % Add the current value to the cumulative sum
            if cumSum > threshold
                times(i:end) = times(i:end) - times(i);
                cumSum = 0;
            end
        end
        shift_cutoff_inds = find(times == 0);

        for k = 1:length(shift_cutoff_inds)
            % if at last index
            if k == length(shift_cutoff_inds)
                shift_bounds = shift_cutoff_inds(k):height(tbl_provider);
            else
                shift_bounds = shift_cutoff_inds(k):shift_cutoff_inds(k+1)-1;
            end
            tbl_shift = tbl_provider(shift_bounds, :);
            tmp_order_in_shift(shift_bounds) = 1:height(tbl_shift);
            % if height(tbl_shift) >= 1
                tmp_time_in_shift(shift_bounds) = hours(tbl_shift.Order_Time - tbl_shift.Order_Time(1));
            % end
            tmp_max_time = hours(tbl_shift.Order_Time(end) - tbl_shift.Order_Time(1));
            tmp_max_time_in_shift(shift_bounds) = repmat(tmp_max_time, numel(shift_bounds), 1);

            % first order
            shift_first_order = tbl_shift.Patient_First_Order;
            if any(shift_first_order) % at least 1 first order for a new patient
                idx_ones = find(shift_first_order); % Find the indices of the '1's in the shift_first_order vector
                % FORWARDS
                b = NaN(size(shift_first_order)); % Calculate the cumulative distance from each element to the nearest '1'               
                for i = 1:length(idx_ones)-1 % Fill the vector b with distances
                    b(idx_ones(i):idx_ones(i+1)-1) = 0:(idx_ones(i+1)-idx_ones(i)-1);
                end
                b(idx_ones(end):end) = 0:(length(shift_first_order)-idx_ones(end));
    
                tmp_order_relative_to_first_forwards(shift_bounds) = b;
    
                % BACKWARDS
                b = NaN(size(shift_first_order));
                for i = 1:length(idx_ones)-1
                    b(idx_ones(i+1):-1:idx_ones(i)+1) = 0:-1:(idx_ones(i)-idx_ones(i+1)+1);
                end
                b(idx_ones(1):-1:1) = 0:-1:-(idx_ones(1)-1);
    
                tmp_order_relative_to_first_backwards(shift_bounds) = b;
            end
    
            % last order
            shift_last_order = tbl_shift.Patient_Last_Order;
            if any(shift_last_order) % at least 1 last order for a new patient
                idx_ones = find(shift_last_order);
                % FORWARDS
                b = NaN(size(shift_last_order)); % Calculate the cumulative distance from each element to the nearest '1'               
                for i = 1:length(idx_ones)-1 % Fill the vector b with distances
                    b(idx_ones(i):idx_ones(i+1)-1) = 0:(idx_ones(i+1)-idx_ones(i)-1);
                end
                b(idx_ones(end):end) = 0:(length(shift_last_order)-idx_ones(end));
    
                tmp_order_relative_to_last_forwards(shift_bounds) = b;
    
                % BACKWARDS
                b = NaN(size(shift_last_order));
                for i = 1:length(idx_ones)-1
                    b(idx_ones(i+1):-1:idx_ones(i)+1) = 0:-1:(idx_ones(i)-idx_ones(i+1)+1);
                end
                b(idx_ones(1):-1:1) = 0:-1:-(idx_ones(1)-1);
    
                tmp_order_relative_to_last_backwards(shift_bounds) = b;
            end
        end
    
        order_relative_to_first_backwards(curr_provider_indices) = tmp_order_relative_to_first_backwards;
        order_relative_to_first_forwards(curr_provider_indices) = tmp_order_relative_to_first_forwards;
        order_relative_to_last_backwards(curr_provider_indices) = tmp_order_relative_to_last_backwards;
        order_relative_to_last_forwards(curr_provider_indices) = tmp_order_relative_to_last_forwards;
        order_in_shift(curr_provider_indices) = tmp_order_in_shift;
        time_in_shift(curr_provider_indices) = tmp_time_in_shift;
        max_time_in_shift(curr_provider_indices) = tmp_max_time_in_shift;
    
    end
end
%%
data_table.("Order_Rel_First_BACKWARDS") = order_relative_to_first_backwards;
data_table.("Order_Rel_First_FORWARDS") = order_relative_to_first_forwards;
data_table.("Order_Rel_Last_BACKWARDS") = order_relative_to_last_backwards;
data_table.("Order_Rel_Last_FORWARDS") = order_relative_to_last_forwards;
data_table.("Order_In_Shift") = order_in_shift;
data_table.("Time_In_Shift") = time_in_shift;
data_table.("Max_Time_In_Shift") = max_time_in_shift;

%%
save('/Users/bib002/Dropbox/research/data/ed_perseveration/data_table_all_order_number.mat','data_table')