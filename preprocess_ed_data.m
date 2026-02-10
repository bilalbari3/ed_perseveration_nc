clear; close all; clc
load('/Users/bib002/Dropbox/research/data/ed_perseveration/tab.mat')
fprintf('Loaded data\n')

% HEADERS
% PatientEncounterID: unique patient identifier
% PrescribingOrAuthorizingProviderID: unique provider identifier
% PrescribingOrAuthorizingSpecialtyDSC: provider specialty
% InstantiatedTimeDTS: time orders released
% ProcedureID: unique order ID
% OrderTypeDSC: category of order
% ProcedureDSC: exact order
% ActivePeopleCountNBR: number of people in ED
% AgeNBR: age
% SexDSC: sex
% RaceDSC: race
% EthnicGroupDSC: ethnicity
% LanguageDSC: language
% FinancialClassDSC: insurance
% EffectiveEpicDepartmentDSC: unclear

%% analyze data
% do this year-by-year to speed it up

% clean up numbers so they're easier to work with
[~,~,tab.PatientEncounterID] = unique(tab.PatientEncounterID);
[~,~,tab.PrescribingOrAuthorizingProviderID] = unique(tab.PrescribingOrAuthorizingProviderID);
[~,~,tab.ProcedureID] = unique(tab.ProcedureID);
[~,~,tab.OrderTypeDSC] = unique(tab.OrderTypeDSC);

num_unique_orders = length(unique(tab.ProcedureID));
num_unique_order_categories = length(unique(tab.OrderTypeDSC));

tab = tab(tab.InstantiatedTimeDTS >= datetime('2023-12-31') & ...
          tab.InstantiatedTimeDTS <= datetime('2025-01-01'), :);
%% batch orders into a vector of orders per unique order time

VariableNames = {'Patient','Provider','Time','Order','Order_Category','Census', ...
                 'Provider_Speciality', ...
                 'Age', ...
                 'Sex', ...
                 'Race', ...
                 'Ethnicity', ...
                 'Language', ...
                 'Insurance', ...
                 'Epic_Department'};
tbl_batch = table(double([]), double([]), NaT(0), {}, {}, double([]), ...
                  categorical([]), ...
                  double([]), categorical([]), categorical([]), categorical([]), ...
                  categorical([]), categorical([]), categorical([]), ...
                  'VariableNames', VariableNames);

i = 1; % go through each entry
while i < height(tab)
    fprintf('On %i of %i\n', i, height(tab))

    patient = tab.PatientEncounterID(i);
    provider = tab.PrescribingOrAuthorizingProviderID(i);
    time = tab.InstantiatedTimeDTS(i);
    census = tab.ActivePeopleCountNBR(i);
    provider_specialty = tab.PrescribingOrAuthorizingSpecialtyDSC(i);
    age = tab.AgeNBR(i);
    sex = tab.SexDSC(i);
    race = tab.RaceDSC(i);
    ethnicity = tab.EthnicGroupDSC(i);
    language = tab.LanguageDSC(i);
    finance = tab.FinancialClassDSC(i);
    department = tab.EffectiveEpicDepartmentDSC(i);
    
    % generate order vector
    order_mask = tab.PatientEncounterID == patient & ...
                       tab.PrescribingOrAuthorizingProviderID == provider & ...
                       tab.InstantiatedTimeDTS == time; % pull out all orders from the batch
    order = tab.ProcedureID(order_mask); % what were the orders?
    order_template = zeros(1, num_unique_orders); % initialize; fill this in with orders that happen at a specific time
    order_template(order) = 1; % assign orders to the order template (to be used for perseveration analysis)
    
    % generate order category vector
    order_category = tab.OrderTypeDSC(order_mask);
    order_category_template = zeros(1, num_unique_order_categories);
    order_category_template(order_category) = 1;

    tbl_tmp = table(patient, provider, time, {order_template}, {order_category_template}, census, ...
                    provider_specialty, ...
                    age, sex, race, ethnicity, language, finance, department, ...
                    'VariableNames', VariableNames);
    tbl_batch = [tbl_batch; tbl_tmp]; % build this up

    i = find(order_mask, 1, 'last') + 1; % iterate to the next line
end

tbl_batch_2024 = tbl_batch;
save('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_2024.mat','tbl_batch_2024','-v7.3')

%% now batch all the data together
clear; close all; clc
load('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_2019.mat')
load('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_2020.mat')
load('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_2021.mat')
load('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_2022.mat')
load('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_2023.mat')
load('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_2024.mat')
tbl_batch = [tbl_batch_2019; tbl_batch_2020; tbl_batch_2021; tbl_batch_2022; tbl_batch_2023; tbl_batch_2024];

save('/Users/bib002/Dropbox/research/data/ed_perseveration/tbl_batch_all.mat','tbl_batch','-v7.3')