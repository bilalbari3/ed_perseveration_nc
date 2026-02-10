load('/Users/bib002/Dropbox/research/data/ed_perseveration/tab.mat')

% clean up numbers so they're easier to work with
[~,~,tab.PatientEncounterID] = unique(tab.PatientEncounterID);
[~,~,tab.PrescribingOrAuthorizingProviderID] = unique(tab.PrescribingOrAuthorizingProviderID);
[~,~,tab.ProcedureID] = unique(tab.ProcedureID);
[~,~,tab.OrderTypeDSC] = unique(tab.OrderTypeDSC);

unique_orders = unique(tab.ProcedureID)';

order_to_cat_LUT = NaN(1, length(unique_orders));
for curr_order = unique_orders
    order_ind = find(tab.ProcedureID == curr_order, 1, 'first');
    order_to_cat_LUT(curr_order) = tab.OrderTypeDSC(order_ind);
end

save('/Users/bib002/Dropbox/research/data/ed_perseveration/order_to_cat_LUT.mat','order_to_cat_LUT')