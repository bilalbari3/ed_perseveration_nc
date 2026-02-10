load('/Users/bib002/Dropbox/research/data/ed_perseveration/tab.mat')
%%
n_oop = []; % orders per patient

for patient = unique(tab.PatientEncounterID')
    tmp_n = sum(tab.PatientEncounterID == patient);
    n_oop = [n_oop; tmp_n];
end

%%
save('/Users/bib002/Dropbox/research/data/ed_perseveration/orders_per_patient.mat', 'n_oop')