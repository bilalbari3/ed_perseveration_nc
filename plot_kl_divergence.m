%% load data
clear; close all; clc
load('/Users/bib002/Dropbox/research/data/ed_perseveration/KL_div_data_table.mat') % dt and dt_struct
saveFigLoc = '/Users/bib002/Documents/git-repositories/ed_perseveration/figures';
myColors = importColors_bb();

%% calculate information-theoretic quantities
bin_edges = 40:1:210;

census1 = dt.Census;
census1(dt.Number_Ordering_Providers ~= 1) = NaN; % remove instances where more than 1 provider placed orders
census2 = dt.Census;
census2(dt.Number_Ordering_Providers <= 1) = NaN;

census1_plot = [];
census2_plot = [];
MI1 = [];
MI2 = [];
HSA1 = [];
HSA2 = [];
KLdiv1 = [];
KLdiv2 = [];
for bin_index = 1:length(bin_edges) - 1
    bin_left = bin_edges(bin_index);
    bin_right = bin_edges(bin_index + 1);
    
    mask_census1 = census1 >= bin_left & census1 < bin_right;
    order1_matrix = dt.Order_Category2(mask_census1, :);

    mask_census2 = census2 >= bin_left & census2 < bin_right;
    order2_matrix = dt.Order_Category2(mask_census2, :);

    Po1 = mean(order1_matrix, "omitnan"); Po1 = Po1./sum(Po1);
    Po2 = mean(order2_matrix, "omitnan"); Po2 = Po2./sum(Po2);

    tmp_cond_ent1 = [];
    tmp_KL_div1 = [];
    for n = 1:size(order1_matrix, 1)
        tmp_orders = order1_matrix(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders); % convert to probabilities
        tmp_cond_ent1(n) = -1*sum(tmp_orders.*log2(tmp_orders), "omitnan");
        tmp_KL_div1(n) = kldiv(1:2, ...
                               tmp_orders + eps, ...
                               Po1 + eps);
    end

    tmp_cond_ent2 = [];
    tmp_KL_div2 = [];
    for n = 1:size(order2_matrix, 1)
        tmp_orders = order2_matrix(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders); % convert to probabilities
        tmp_cond_ent2(n) = -1*sum(tmp_orders.*log2(tmp_orders), "omitnan");
        tmp_KL_div2(n) = kldiv(1:2, ...
                               tmp_orders + eps, ...
                               Po2 + eps);
    end

    census1_plot = [census1_plot; mean([bin_left bin_right])];
    census2_plot = [census2_plot; mean([bin_left bin_right])];

    H_Po1 = -1*sum(Po1.*log2(Po1),'omitnan');
    MI1 = [MI1; H_Po1 - mean(tmp_cond_ent1)];
    HSA1 = [HSA1; mean(tmp_cond_ent1)];
    KLdiv1 = [KLdiv1; mean(tmp_KL_div1)];
    
    H_Po2 = -1*sum(Po2.*log2(Po2),'omitnan');
    MI2 = [MI2; H_Po2 - mean(tmp_cond_ent2)];
    HSA2 = [HSA2; mean(tmp_cond_ent2)];
    KLdiv2 = [KLdiv2; mean(tmp_KL_div2)];
end

%% plot info-theoretic quantities - 1 provider only
n_bins = 11;
plot_color1 = myColors.reddishPurple;

f_info = figure('Position',[411    55   400   845]); 
h_MI = subplot(311); hold on; 
ylabel(sprintf('Mutual information\nI(orders ; patients)'))
h_HSA = subplot(312); hold on; 
ylabel('H(orders | patient)')
h_KLdiv = subplot(313); hold on; 
ylabel(sprintf('KL divergence\nD_{KL}(P(orders | patients) || P(orders))')); xlabel('Census')

t(1) = plotBinned(census1_plot, MI1, n_bins, h_MI, 'Color', plot_color1);
ylim([0.08 0.11])

t(2) = plotBinned(census1_plot, HSA1, n_bins, h_HSA, 'Color', plot_color1);
ylim([0.88 0.91])

t(3) = plotBinned(census1_plot, KLdiv1, n_bins, h_KLdiv, 'Color', plot_color1);
ylim([0.08 0.11])


set([h_MI h_HSA h_KLdiv], 'fontsize', 18, 'tickdir', 'out')

saveFigureIteration(f_info, '/Users/bib002/Documents/git-repositories/ed_perseveration/figures', ...
                    'info_measures_1provider')

%% info-theoretic quantities (split by # providers)
n_bins = 11;
plot_color1 = myColors.reddishPurple;
plot_color2 = myColors.bluishGreen;

f_info = figure('Position',[179    21   784   824]); 
h_MI = subplot(321); hold on; 
ylabel(sprintf('Mutual information\nI(orders ; patients)'))
h_HSA = subplot(323); hold on; 
ylabel('H(orders | patient)')
h_KLdiv = subplot(325); hold on; 
ylabel(sprintf('KL divergence\nD_{KL}(P(orders | patients) || P(orders))')); xlabel('Census')

t(1) = plotBinned(census1_plot, MI1, n_bins, h_MI, 'Color', plot_color1);
t(2) = plotBinned(census2_plot, MI2, n_bins, h_MI, 'Color', plot_color2);
legend(t,{'1 Provider','>1 Provider'})

plotBinned(census1_plot, HSA1, n_bins, h_HSA, 'Color', plot_color1);
plotBinned(census2_plot, HSA2, n_bins, h_HSA, 'Color', plot_color2);

plotBinned(census1_plot, KLdiv1, n_bins, h_KLdiv, 'Color', plot_color1);
plotBinned(census2_plot, KLdiv2, n_bins, h_KLdiv, 'Color', plot_color2);

h_MI_hist = subplot(322); hold on
bin_edges = linspace(0,0.3,20);
histogram(MI1, bin_edges, 'Normalization', 'probability', 'FaceColor', plot_color1, 'EdgeColor', 'none')
histogram(MI2, bin_edges, 'Normalization', 'probability', 'FaceColor', plot_color2, 'EdgeColor', 'none')
xlabel('I(orders ; patients)')
xlim([0.05 0.25])

h_HSA_hist = subplot(324); hold on
bin_edges = linspace(0.7, 1, 20);
histogram(HSA1, bin_edges, 'Normalization', 'probability', 'FaceColor', plot_color1, 'EdgeColor', 'none')
histogram(HSA2, bin_edges, 'Normalization', 'probability', 'FaceColor', plot_color2, 'EdgeColor', 'none')
xlabel('Conditional entropy')
xlim([0.75 0.95])

h_KLdiv_hist = subplot(326); hold on
bin_edges = linspace(0,0.3,20);
histogram(KLdiv1, bin_edges, 'Normalization', 'probability', 'FaceColor', plot_color1, 'EdgeColor', 'none')
histogram(KLdiv2, bin_edges, 'Normalization', 'probability', 'FaceColor', plot_color2, 'EdgeColor', 'none')
xlabel('KL divergence')
xlim([0.05 0.25])


set([h_MI h_HSA h_KLdiv h_MI_hist h_HSA_hist h_KLdiv_hist], 'fontsize', 18,'tickdir','out')
ylabel([h_MI_hist h_HSA_hist h_KLdiv_hist], 'Probability')
saveFigureIteration(f_info, '/Users/bib002/Documents/git-repositories/ed_perseveration/figures', ...
                    'info_measures_1plusprovider')

%% fit regressions
dt_regression = table(normalize(MI1), normalize(HSA1), normalize(KLdiv1), normalize(census1_plot), ...
                      normalize(MI2), normalize(HSA2), normalize(KLdiv2), normalize(census2_plot), ...
                     'VariableNames', ...
                    {'MI1','HSA1','KLdiv1','Census1', ...
                     'MI2','HSA2','KLdiv2','Census2'});
mod_MI1 = fitlm(dt_regression, 'MI1 ~ Census1');
mod_HSA1 = fitlm(dt_regression, 'HSA1 ~ Census1');
mod_KLdiv1 = fitlm(dt_regression, 'KLdiv1 ~ Census1');

%% oxloc data to test calibration of priors
rng(1)
load('/Users/bib002/Dropbox/research/data/ed_perseveration/oxloc.mat'); % load oxloc data
oxloc_MA = oxloc(strcmp(oxloc.RegionCode, "US_MA"), :);

% pick data
data_source = 'COVID_death'; % COVID_cases, COVID_death, COVID_CHI
switch data_source
    case 'COVID_cases'
        covid_cases = [0; diff(oxloc_MA.ConfirmedCases)];
        covid_cases = smooth(covid_cases, 50, 'sgolay', 2); % smooth cases/deaths
        covid_cases(covid_cases < 0) = 0;
    case 'COVID_death'
        covid_cases = [0; diff(oxloc_MA.ConfirmedDeaths)]; covid_cases([182 232]) = 0; % highly negative in deaths
        covid_cases = smooth(covid_cases, 50, 'sgolay', 2); % smooth cases/deaths
        covid_cases(covid_cases < 0) = 0;
    case 'COVID_CHI'
        covid_cases = oxloc_MA.ContainmentHealthIndex_WeightedAverage;
end

% compare 1st and 3rd tercile
covid_median = median(covid_cases);
covid_thresh = prctile(covid_cases, [33 66]);
covid_thresh_low = covid_thresh(1);
covid_thresh_high = covid_thresh(2);

f_COVID = figure; hold on
plot(oxloc_MA.Date, covid_cases, 'linewidth', 2, 'Color', [0 0 0])
mask_covid_low = covid_cases <= covid_thresh_low;
mask_covid_high = covid_cases > covid_thresh_high;
plot(oxloc_MA.Date(mask_covid_high), covid_cases(mask_covid_high), 'ro')
plot(oxloc_MA.Date(mask_covid_low), covid_cases(mask_covid_low), 'bo')

set(gca,'fontsize',18,'tickdir','out')

saveFigureIteration(f_COVID, '/Users/bib002/Documents/git-repositories/ed_perseveration/figures', ...
                    strcat('timeseries_', data_source))

%%
% shift dates to start of day for comparison
dt_dates = dateshift(dt.Time, 'start', 'day');

% shuffle the data; 1st index is unshuffled; all others are shuffled
KLdiv_diff = [];
KLdiv_COVID_high_cell = [];
KLdiv_COVID_low_cell = [];
for i = 1:101
    fprintf('On %i\n', i)

    % construct masks for high and low COVID dates
    mask_covid_low = covid_cases <= covid_thresh_low;
    mask_covid_high = covid_cases > covid_thresh_high;
    
    % shuffle the data
    if i ~= 1
        mask_covid_low = mask_covid_low(randperm(length(mask_covid_low))); % shuffle the data
        mask_covid_high = mask_covid_high(randperm(length(mask_covid_high))); % shuffle the data
    else
        fprintf('Not shuffling\n')
    end
    
    % mask for low covid dates
    covid_low_dates = (oxloc_MA.Date(mask_covid_low));
    m_covid_low_dates = false(length(dt_dates), 1);
    for curr_date = covid_low_dates'
        m_covid_low_dates(dt_dates == curr_date) = true;
    end
    
    % mask for high covid dates
    covid_high_dates = (oxloc_MA.Date(mask_covid_high));
    m_covid_high_dates = false(length(dt_dates), 1);
    for curr_date = covid_high_dates'
        m_covid_high_dates(dt_dates == curr_date) = true;
    end
    % only look at dates w/ 1 ordering provider
    m_covid_low_dates(dt.Number_Ordering_Providers ~= 1) = false;
    m_covid_high_dates(dt.Number_Ordering_Providers ~= 1) = false;
    
    orders_COVID_low = dt.Order_Category2(m_covid_low_dates, :);
    orders_COVID_high = dt.Order_Category2(m_covid_high_dates, :);
    
    Po_COVID_low = mean(orders_COVID_low,"omitnan"); Po_COVID_low = Po_COVID_low./sum(Po_COVID_low);
    Po_COVID_high = mean(orders_COVID_high,"omitnan"); Po_COVID_high = Po_COVID_high./sum(Po_COVID_high);
    
    KLdiv_COVID_low = [];
    KLdiv_COVID_high = [];
    
    for n = 1:size(orders_COVID_low, 1)
        tmp_orders = orders_COVID_low(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders);
        KLdiv_COVID_low(n) = kldiv(1:2, tmp_orders + eps, Po_COVID_high + eps);
    end
    for n = 1:size(orders_COVID_high, 1)
        tmp_orders = orders_COVID_high(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders);
        KLdiv_COVID_high(n) = kldiv(1:2, tmp_orders + eps, Po_COVID_high + eps);
    end

    KLdiv_diff = [KLdiv_diff; mean(KLdiv_COVID_high) - mean(KLdiv_COVID_low)];
    KLdiv_COVID_high_cell{i} = KLdiv_COVID_high;
    KLdiv_COVID_low_cell{i} = KLdiv_COVID_low;
end
%% save data
data_file = strcat('oxloc_', data_source, '.mat');
save(fullfile('/Users/bib002/Dropbox/research/data/ed_perseveration', data_file), ...
    'KLdiv_diff','KLdiv_COVID_high_cell','KLdiv_COVID_low_cell')
fprintf('Saved %s\n', data_file)
%% load data
data_file = 'COVID_death'; % COVID_cases, COVID_death, COVID_CHI
load(fullfile('/Users/bib002/Dropbox/research/data/ed_perseveration/', ...
     strcat('oxloc_', data_file, '.mat')))
%%
f_COVID = figure('Position',[161     1   507   865]); 
h_KLdiv = subplot(211); hold on
errorbar(1, mean(KLdiv_COVID_low_cell{1}), sem(KLdiv_COVID_low_cell{1}), '.k', 'linewidth', 2,'MarkerSize',30)
errorbar(2, mean(KLdiv_COVID_high_cell{1}), sem(KLdiv_COVID_high_cell{1}), '.k', 'linewidth', 2,'MarkerSize',30)
xlim([0 3])
% switch data_file
%     case 'COVID_death'
%         ylim([0.084 0.094])
%     case 'COVID_CHI'
%         ylim([0.084 0.094])
% end
set(gca,'xtick',1:2,'xticklabel',{'COVID low', 'COVID high'},'fontsize',18,'tickdir','out')
xlabel('Condition')
ylabel(sprintf('KL Div\nP(order | patient) || P(order | COVID low)'))

h_shuffle = subplot(212); hold on
histogram(KLdiv_diff(2:end), 'normalization', 'probability', 'EdgeColor', 'none')
plot([KLdiv_diff(1) KLdiv_diff(1)], [0 0.25], 'k', 'linewidth', 2)
set(gca,'fontsize',18,'tickdir','out')
xlabel('KL Div')
ylabel('Probability')
legend({'Shuffle','Data'},'location','best')
% switch data_file
%     case 'COVID_death'
%         xlim([-3 7]*10^-3)
%     case 'COVID_CHI'
%         xlim([-3 7]*10^-3)
% end

saveFigureIteration(f_COVID, '/Users/bib002/Documents/git-repositories/ed_perseveration/figures', ...
    data_file)

%% compare high vs low for both priors
rng(1)
load('/Users/bib002/Dropbox/research/data/ed_perseveration/oxloc.mat'); % load oxloc data
oxloc_MA = oxloc(strcmp(oxloc.RegionCode, "US_MA"), :);

% pick data
data_source = 'COVID_cases'; % COVID_cases, COVID_death, COVID_CHI
switch data_source
    case 'COVID_cases'
        covid_cases = [0; diff(oxloc_MA.ConfirmedCases)];
        covid_cases = smooth(covid_cases, 50, 'sgolay', 2); % smooth cases/deaths
        covid_cases(covid_cases < 0) = 0;
    case 'COVID_death'
        covid_cases = [0; diff(oxloc_MA.ConfirmedDeaths)]; covid_cases([182 232]) = 0; % highly negative in deaths
        covid_cases = smooth(covid_cases, 50, 'sgolay', 2); % smooth cases/deaths
        covid_cases(covid_cases < 0) = 0;
    case 'COVID_CHI'
        covid_cases = oxloc_MA.ContainmentHealthIndex_WeightedAverage;
end

% calculate thresholds
covid_thresh = prctile(covid_cases, [33 66]);
covid_thresh_low = covid_thresh(1);
covid_thresh_high = covid_thresh(2);

% shift dates to start of day for comparison
dt_dates = dateshift(dt.Time, 'start', 'day');

% shuffle the data; 1st index is unshuffled; all others are shuffled
KLdiv_diff_low = []; % prior high minus prior low
KLdiv_COVID_low_priorhigh_cell = [];
KLdiv_COVID_low_priorlow_cell = [];

KLdiv_diff_high = []; % prior high minus prior low
KLdiv_COVID_high_priorhigh_cell = [];
KLdiv_COVID_high_priorlow_cell = [];
for i = 1:1001
    fprintf('On %i\n', i)

    % construct masks for high and low COVID dates
    mask_covid_low = covid_cases <= covid_thresh_low;
    mask_covid_high = covid_cases > covid_thresh_high;
    
    % shuffle the data
    if i ~= 1
        mask_covid_low = mask_covid_low(randperm(length(mask_covid_low))); % shuffle the data
        mask_covid_high = mask_covid_high(randperm(length(mask_covid_high))); % shuffle the data
    else
        fprintf('Not shuffling\n')
    end
    
    % mask for low covid dates
    covid_low_dates = (oxloc_MA.Date(mask_covid_low));
    m_covid_low_dates = false(length(dt_dates), 1);
    for curr_date = covid_low_dates'
        m_covid_low_dates(dt_dates == curr_date) = true;
    end
    
    % mask for high covid dates
    covid_high_dates = (oxloc_MA.Date(mask_covid_high));
    m_covid_high_dates = false(length(dt_dates), 1);
    for curr_date = covid_high_dates'
        m_covid_high_dates(dt_dates == curr_date) = true;
    end
    % only look at dates w/ 1 ordering provider
    m_covid_low_dates(dt.Number_Ordering_Providers ~= 1) = false;
    m_covid_high_dates(dt.Number_Ordering_Providers ~= 1) = false;
    
    orders_COVID_low = dt.Order_Category2(m_covid_low_dates, :);
    orders_COVID_high = dt.Order_Category2(m_covid_high_dates, :);
    
    Po_COVID_low = mean(orders_COVID_low,"omitnan"); Po_COVID_low = Po_COVID_low./sum(Po_COVID_low);
    Po_COVID_high = mean(orders_COVID_high,"omitnan"); Po_COVID_high = Po_COVID_high./sum(Po_COVID_high);
    
    % COVID low dates first
    KLdiv_COVID_low_priorlow = [];
    KLdiv_COVID_low_priorhigh = [];
    for n = 1:size(orders_COVID_low, 1)
        tmp_orders = orders_COVID_low(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders);
        KLdiv_COVID_low_priorlow(n) = kldiv(1:2, tmp_orders + eps, Po_COVID_low + eps);
    end
    for n = 1:size(orders_COVID_low, 1)
        tmp_orders = orders_COVID_low(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders);
        KLdiv_COVID_low_priorhigh(n) = kldiv(1:2, tmp_orders + eps, Po_COVID_high + eps);
    end    
    
    KLdiv_diff_low = [KLdiv_diff_low; mean(KLdiv_COVID_low_priorhigh) - mean(KLdiv_COVID_low_priorlow)];
    KLdiv_COVID_low_priorhigh_cell{i} = KLdiv_COVID_low_priorhigh;
    KLdiv_COVID_low_priorlow_cell{i} = KLdiv_COVID_low_priorlow;

    % COVID high dates next
    KLdiv_COVID_high_priorlow = [];
    KLdiv_COVID_high_priorhigh = [];
    for n = 1:size(orders_COVID_high, 1)
        tmp_orders = orders_COVID_high(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders);
        KLdiv_COVID_high_priorlow(n) = kldiv(1:2, tmp_orders + eps, Po_COVID_low + eps);
    end
    for n = 1:size(orders_COVID_high, 1)
        tmp_orders = orders_COVID_high(n, :);
        tmp_orders = tmp_orders./sum(tmp_orders);
        KLdiv_COVID_high_priorhigh(n) = kldiv(1:2, tmp_orders + eps, Po_COVID_high + eps);
    end    

    KLdiv_diff_high = [KLdiv_diff_high; mean(KLdiv_COVID_high_priorhigh) - mean(KLdiv_COVID_high_priorlow)];
    KLdiv_COVID_high_priorhigh_cell{i} = KLdiv_COVID_high_priorhigh;
    KLdiv_COVID_high_priorlow_cell{i} = KLdiv_COVID_high_priorlow;
end

%%
data_file = strcat('oxloc_', data_source, '_n1000.mat');
save(fullfile('/Users/bib002/Dropbox/research/data/ed_perseveration', data_file), ...
    'KLdiv_diff_low','KLdiv_COVID_low_priorhigh_cell','KLdiv_COVID_low_priorlow_cell', ...
    'KLdiv_diff_high','KLdiv_COVID_high_priorhigh_cell','KLdiv_COVID_high_priorlow_cell')
fprintf('Saved %s\n', data_file)
%% load data
load('/Users/bib002/Dropbox/research/data/ed_perseveration/oxloc_COVID_cases_n1000.mat')
%% plot
analysis = "high";
switch analysis
    case "high"
        KLdiv_diff = KLdiv_diff_high;
        KLdiv_COVID_low_cell = KLdiv_COVID_high_priorlow_cell;
        KLdiv_COVID_high_cell = KLdiv_COVID_high_priorhigh_cell;
    case "low"
        KLdiv_diff = KLdiv_diff_low;
        KLdiv_COVID_low_cell = KLdiv_COVID_low_priorlow_cell;
        KLdiv_COVID_high_cell = KLdiv_COVID_low_priorhigh_cell;
end

dat = [KLdiv_COVID_low_cell{1}' KLdiv_COVID_high_cell{1}'];
[se, m] = wse(dat, 2);

f_COVID = figure('Position',[161     1   313   734]); 
h_KLdiv = subplot(211); hold on
errorbar(1, m(1), se(1), '.k', 'linewidth', 2,'MarkerSize',30)
errorbar(2, m(2), se(2), '.k', 'linewidth', 2,'MarkerSize',30)
xlim([0.5 2.5])
switch analysis
    case "low"
        ylim([0.0890 0.0898])
    case "high"
        ylim([0.0884 0.0892])
end

set(gca,'xtick',1:2,'xticklabel',{'COVID low', 'COVID high'},'fontsize',18,'tickdir','out')
xlabel('Condition')
ylabel(sprintf('KL Div\nP(order | patient) || P(order | COVID low)'))

h_shuffle = subplot(212); hold on
histogram(KLdiv_diff(2:end), 10, 'normalization', 'probability', 'EdgeColor', 'none')
plot([KLdiv_diff(1) KLdiv_diff(1)], [0 0.25], 'k', 'linewidth', 2)
set(gca,'fontsize',18,'tickdir','out')
xlabel('KL Div')
ylabel('Probability')
legend({'Shuffle','Data'},'location','best')
switch analysis
    case "low"
        xlim([0 6e-4])
        ylim([0 0.4])
    case "high"
        xlim([-6e-4 0])
end

saveFigureIteration(f_COVID, '/Users/bib002/Documents/git-repositories/ed_perseveration/figures', ...
    char(strcat("COVID_cases_", analysis)))

%%
d1 = KLdiv_COVID_low_priorlow_cell{1};
d2 = KLdiv_COVID_high_priorlow_cell{1};
% dat = [KLdiv_COVID_low_priorlow_cell{1}' KLdiv_COVID_high_priorlow_cell{1}'];
% [se, m] = wse(dat, 2);

f_COVID = figure('Position',[161     1   507   865]); 
h_KLdiv = subplot(211); hold on
errorbar(1, mean(d1), sem(d1), '.k', 'linewidth', 2,'MarkerSize',30)
errorbar(2, mean(d2), sem(d2), '.k', 'linewidth', 2,'MarkerSize',30)
xlim([0.5 2.5])