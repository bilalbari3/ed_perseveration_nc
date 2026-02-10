%% 2 set size
rng(3) % 3 works well
n_trials = 1e3;
Ps = ones(1, 2)/2;
b = logspace(-1,3,1e2);

R_2 = NaN(n_trials, length(b));
V_2 = R_2;
HSA_2 = R_2;
KLdiv_2 = R_2;
for i = 1:n_trials
    escape_flag = false;
    if mod(i,100) == 0
        fprintf('Set size 2: On %i\n', i)
    end
    while escape_flag == false % ensure optimal marginal action distribution isn't [1 0]
        Q = rand(2);
        [~, Qmax_index] = max(Q, [], 2);
        if isscalar(unique(Qmax_index)) == false
            escape_flag = true;
        end
    end

    [R,V,Pa,Psa] = blahut_arimoto(Ps,Q,b);
    R_2(i,:) = R;
    V_2(i,:) = V;

    H_Pa = -1*sum(Pa'.*log2(Pa'));
    tmp_HSA = [];
    tmp_KLdiv = [];
    for j = 1:length(Psa)
        tmp_HSA(j) = mean(-1*sum(Psa{j}'.*log2(Psa{j}')));
        tmp_KLdiv(j) = mean([kldiv(1:2, Psa{j}(1,:) + eps, Pa(j,:) + eps) ...
                             kldiv(1:2, Psa{j}(2,:) + eps, Pa(j,:) + eps)]);
    end
    HSA_2(i,:) = tmp_HSA;
    KLdiv_2(i,:) = tmp_KLdiv;

end

% 4 set size
Ps = ones(1, 4)/4;
R_4 = NaN(n_trials, length(b));
V_4 = R_4;
HSA_4 = R_4;
KLdiv_4 = R_4;
for i = 1:n_trials
    escape_flag = false;
    if mod(i,100) == 0
        fprintf('Set size 4: On %i\n', i)
    end
    while escape_flag == false
        Q = rand(4,2);
        [~, Qmax_index] = max(Q, [], 2);
        if isscalar(unique(Qmax_index)) == false
            escape_flag = true;
        end
    end
    [R,V,Pa,Psa] = blahut_arimoto(Ps,Q,b);

    R_4(i,:) = R;
    V_4(i,:) = V;

    H_Pa = -1*sum(Pa'.*log2(Pa'));
    tmp_HSA = [];
    tmp_KLdiv = [];
    for j = 1:length(Psa)
        tmp_HSA(j) = mean(-1*sum(Psa{j}'.*log2(Psa{j}')));
        tmp_KLdiv(j) = mean([kldiv(1:2, Psa{j}(1,:) + eps, Pa(j,:) + eps) ...
                             kldiv(1:2, Psa{j}(2,:) + eps, Pa(j,:) + eps) ...
                             kldiv(1:2, Psa{j}(3,:) + eps, Pa(j,:) + eps) ...
                             kldiv(1:2, Psa{j}(4,:) + eps, Pa(j,:) + eps)]);
    end
    HSA_4(i,:) = tmp_HSA;
    KLdiv_4(i,:) = tmp_KLdiv;

end
fprintf('Fin\n')

R_2(R_2(:, end) < 1 - 1e-6, :) = [];
R_4(R_4(:, end) < 1 - 1e-6, :) = []; % exclude ones that didn't max out

%% differing policy complexity
I_low = [];
beta_low = [];
ind_I_low = [];
HSA_low = [];
KLdiv_low = [];

I_high = [];
beta_high = [];
ind_I_high = [];
HSA_high = [];
KLdiv_high = [];
for i = 1:size(R_2, 1)
    [~, tmp_ind_I_low] = min(abs(R_2(i,:)-0.16)); % low is 0.16 bits
    I_low = [I_low; R_2(i, tmp_ind_I_low)];
    beta_low = [beta_low; b(tmp_ind_I_low)];
    ind_I_low = [ind_I_low; tmp_ind_I_low];
    HSA_low = [HSA_low; HSA_2(i, tmp_ind_I_low)];
    KLdiv_low = [KLdiv_low; KLdiv_2(i, tmp_ind_I_low)];


    [~, tmp_ind_I_high] = min(abs(R_2(i,:)-0.9966)); % high is 0.99 bits
    I_high = [I_high; R_2(i, tmp_ind_I_high)];
    beta_high = [beta_high; b(tmp_ind_I_high)];
    ind_I_high = [ind_I_high; tmp_ind_I_high];
    HSA_high = [HSA_high; HSA_2(i, tmp_ind_I_high)];
    KLdiv_high = [KLdiv_high; KLdiv_2(i, tmp_ind_I_high)];
end



h_f = figure('Position',[38   231   709   635]); 
% h_RD = subplot(221); hold on
% plot(mean(R_2), normalize(mean(V_2),'range'), 'linewidth', 4, 'Color', myColors.skyBlue_bright)
% xlabel('Policy complexity (bits)'); ylabel('Utility (normalized)')
% set(h_RD,'ytick',0:0.5:1)

h_RB = subplot(221); hold on
plot(mean(R_2), b, 'linewidth', 4, 'Color', myColors.skyBlue_bright)
xlabel('Policy complexity (bits)'); ylabel('Beta')
ylim([0 30])

h_optBeta = subplot(222); hold on
bar(2, mean(beta_low), 'FaceColor',myColors.skyBlue_bright)
errorbar(2, mean(beta_low), sem(beta_low), 'linewidth', 2, 'Color', myColors.black)
bar(1, mean(beta_high), 'FaceColor',myColors.skyBlue_dull)
errorbar(1, mean(beta_high), sem(beta_high), 'linewidth', 2, 'Color', myColors.black)
xlim([0.5 2.5])
xlabel('Policy complexity')
ylabel('Beta')


h_HSA = subplot(223); hold on
plot(2, mean(HSA_low),'.','Color',myColors.skyBlue_bright,'MarkerSize',20)
errorbar(2, mean(HSA_low), sem(HSA_low), 'linewidth', 2, 'Color', myColors.skyBlue_bright)
plot(1, mean(HSA_high),'.','Color',myColors.skyBlue_dull,'MarkerSize',20)
errorbar(1, mean(HSA_high), sem(HSA_high), 'linewidth', 2, 'Color', myColors.skyBlue_dull)
xlim([0.5 2.5])
xlabel('Policy complexity')
ylabel('Conditional entropy')

h_KLdiv = subplot(224); hold on
plot(2, mean(KLdiv_low),'.','Color',myColors.skyBlue_bright,'MarkerSize',20)
errorbar(2, mean(KLdiv_low), sem(KLdiv_low), 'linewidth', 2, 'Color', myColors.skyBlue_bright)
plot(1, mean(KLdiv_high),'.','Color',myColors.skyBlue_dull,'MarkerSize',20)
errorbar(1, mean(KLdiv_high), sem(KLdiv_high), 'linewidth', 2, 'Color', myColors.skyBlue_dull)
xlim([0.5 2.5])
ylim([0.2 0.9])
xlabel('Policy complexity')
ylabel('KL divergence')

set([h_RB h_optBeta h_HSA h_KLdiv],'fontsize',18,'tickdir','out')
set([h_optBeta h_HSA h_KLdiv], 'xtick', 1:2, 'xticklabel', {'High','Low'})

%% cognitive load figure
h_f = figure('Position',[38   231   709   635]); 
h_RD = subplot(221); hold on
plot(mean(R_2), normalize(mean(V_2),'range'), 'linewidth', 4, 'Color', myColors.skyBlue_bright)
plot(mean(R_4), normalize(mean(V_4),'range'), 'linewidth', 4, 'Color', myColors.skyBlue_dull)
% plot(mean(R_2), mean(V_2), 'linewidth', 2)
% plot(mean(R_4), mean(V_4), 'linewidth', 2)
xlabel('Policy complexity (bits)'); ylabel('Utility (normalized)')
set(h_RD,'ytick',0:0.5:1)
legend({'Set size 2','Set size 4'},'location','best')

h_RB = subplot(222); hold on
plot(mean(R_2), b, 'linewidth', 4, 'Color', myColors.skyBlue_bright)
plot(mean(R_4), b, 'linewidth', 4, 'Color', myColors.skyBlue_dull)
xlabel('Policy complexity (bits)'); ylabel('Beta')
ylim([0 30])

h_HSA = subplot(223); hold on
plot(1, mean(HSA_2(:),"omitnan"),'.','Color',myColors.skyBlue_bright,'MarkerSize',20)
errorbar(1, mean(HSA_2(:), "omitnan"), sem(HSA_2(:)), 'linewidth', 2, 'Color', myColors.skyBlue_bright)
plot(2, mean(HSA_4(:),"omitnan"),'.','Color',myColors.skyBlue_dull,'MarkerSize',20)
errorbar(2, mean(HSA_4(:), "omitnan"), sem(HSA_4(:)), 'linewidth', 2, 'Color', myColors.skyBlue_dull)
plot([1 2],[mean(HSA_2(:), "omitnan") mean(HSA_4(:), "omitnan")], 'linewidth', 2, 'Color', myColors.skyBlue)
xlim([0.5 2.5])
xlabel('Cognitive load')
ylabel('Conditional entropy')

h_KLdiv = subplot(224); hold on
plot(1, mean(KLdiv_2(:),"omitnan"),'.','Color',myColors.skyBlue_bright,'MarkerSize',20)
errorbar(1, mean(KLdiv_2(:), "omitnan"), sem(KLdiv_2(:)), 'linewidth', 2, 'Color', myColors.skyBlue_bright)
plot(2, mean(KLdiv_4(:),"omitnan"),'.','Color',myColors.skyBlue_dull,'MarkerSize',20)
errorbar(2, mean(KLdiv_4(:), "omitnan"), sem(KLdiv_4(:)), 'linewidth', 2, 'Color', myColors.skyBlue_dull)
plot([1 2],[mean(KLdiv_2(:), "omitnan") mean(KLdiv_4(:), "omitnan")], 'linewidth', 2, 'Color', myColors.skyBlue)
xlim([0.5 2.5])
xlabel('Cognitive load')
ylabel('KL divergence')

set([h_RD h_RB h_HSA h_KLdiv],'fontsize',18,'tickdir','out')
set([h_HSA h_KLdiv], 'xtick', 1:2, 'xticklabel', {'Low','High'})

saveFigureIteration(h_f, ...
                    '/Users/bib002/Documents/git-repositories/ed_perseveration/figures', ...
                    'theory')

%% cognitive load figure (bar)
h_f = figure('Position',[38   231   709   635]); 
h_RD = subplot(221); hold on
plot(mean(R_2), normalize(mean(V_2),'range'), 'linewidth', 4, 'Color', myColors.skyBlue_bright)
plot(mean(R_4), normalize(mean(V_4),'range'), 'linewidth', 4, 'Color', myColors.skyBlue_dull)
% plot(mean(R_2), mean(V_2), 'linewidth', 2)
% plot(mean(R_4), mean(V_4), 'linewidth', 2)
xlabel('Policy complexity (bits)'); ylabel('Utility (normalized)')
set(h_RD,'ytick',0:0.5:1)
legend({'Set size 2','Set size 4'},'location','best')

h_RB = subplot(222); hold on
plot(mean(R_2), b, 'linewidth', 4, 'Color', myColors.skyBlue_bright)
plot(mean(R_4), b, 'linewidth', 4, 'Color', myColors.skyBlue_dull)
xlabel('Policy complexity (bits)'); ylabel('Beta')
ylim([0 30])

h_HSA = subplot(223); hold on
bar(1, mean(HSA_2(:),"omitnan"),'FaceColor',myColors.skyBlue_bright,'EdgeColor','none')
errorbar(1, mean(HSA_2(:), "omitnan"), sem(HSA_2(:)), 'linewidth', 2, 'Color', myColors.black)
bar(2, mean(HSA_4(:),"omitnan"),'FaceColor',myColors.skyBlue_dull,'EdgeColor','none')
errorbar(2, mean(HSA_4(:), "omitnan"), sem(HSA_4(:)), 'linewidth', 2, 'Color', myColors.black)
xlim([0.5 2.5])
ylim([0.13 0.16])
xlabel('Cognitive load')
ylabel('Conditional entropy')

h_KLdiv = subplot(224); hold on
bar(1, mean(KLdiv_2(:),"omitnan"),'FaceColor',myColors.skyBlue_bright,'EdgeColor','none')
errorbar(1, mean(KLdiv_2(:), "omitnan"), sem(KLdiv_2(:)), 'linewidth', 2, 'Color', myColors.black)
bar(2, mean(KLdiv_4(:),"omitnan"), 'FaceColor',myColors.skyBlue_dull,'EdgeColor','none')
errorbar(2, mean(KLdiv_4(:), "omitnan"), sem(KLdiv_4(:)), 'linewidth', 2, 'Color', myColors.black)
xlim([0.5 2.5])
ylim([0.42 0.48])
xlabel('Cognitive load')
ylabel('KL divergence')

set([h_RD h_RB h_HSA h_KLdiv],'fontsize',18,'tickdir','out')
set([h_HSA h_KLdiv], 'xtick', 1:2, 'xticklabel', {'Low','High'})

saveFigureIteration(h_f, ...
                    '/Users/bib002/Documents/git-repositories/ed_perseveration/figures', ...
                    'theory_bar')

%% choice repetition
rng(2)
% Q = [1 0; % [ patient x order ]
%      0 1];
% Q = [0.0104    0.4958
%      0.5019    0.1338];
n_repeat = 1e2;
all_persev = [];
allQ = [];
for curr_repeat = 1:n_repeat
    if mod(curr_repeat, 100) == 0
        fprintf('On %i of %i\n', curr_repeat, n_repeat)
    end
    
    escape_flag = false;
    while escape_flag == false % ensure optimal marginal action distribution isn't [1 0]
        Q = rand(2,2);
        [~, Qmax_index] = max(Q, [], 2);
        if isscalar(unique(Qmax_index)) == false
            escape_flag = true;
        end
    end
    % Q(2,:) = Q(1,end:-1:1);
    % allQ(:,:,curr_repeat) = Q;
    
    persev = [];
    % beta_set = linspace(0.1,10,10);
    beta_set = [100 1]; % [low load, high load]
    n_trials = 1e3;
    patients = randi(2,1,n_trials);
    alpha_persev = 0.03;
    
    for beta = beta_set
        Pa = [0.5 0.5];
        choices = [];
        for i = 1:n_trials
            Qdiff = Q(patients(i), 2) - Q(patients(i), 1);
            p_choice = exp(beta*Q(patients(i), 2) + log(Pa(2))) ./ ...
                      (exp(beta*Q(patients(i), 1) + log(Pa(1))) + exp(beta*Q(patients(i), 2) + log(Pa(2))));
            choices(i) = binornd(1, p_choice);
            Pa(1) = Pa(1) + alpha_persev*((1 - p_choice) - Pa(1));
            Pa(2) = Pa(2) + alpha_persev*(p_choice - Pa(2));
        end
        persev = [persev mean(diff(choices) == 0)];
    end
    all_persev = [all_persev; persev];
end
mean(all_persev)
%%
f_persev = figure('Position',[456   474   392   185]); hold on
plot(1, mean(all_persev(:, 1),"omitnan"),'.','Color',myColors.skyBlue_bright,'MarkerSize',20)
errorbar(1, mean(all_persev(:, 1), "omitnan"), sem(all_persev(:, 1)), 'linewidth', 4, 'Color', myColors.skyBlue_bright)
plot(2, mean(all_persev(:, 2),"omitnan"),'.','Color',myColors.skyBlue_dull,'MarkerSize',20)
errorbar(2, mean(all_persev(:, 2), "omitnan"), sem(all_persev(:, 1)), 'linewidth', 4, 'Color', myColors.skyBlue_dull)
plot([1 2],mean(all_persev), 'linewidth', 2, 'Color', myColors.skyBlue)
xlim([0.5 2.5])
% ylim([0.13 0.16])
xlabel('Cognitive load')
ylabel('P(repeat order)')

set(gca,'fontsize',18,'tickdir','out')
set(gca, 'xtick', 1:2, 'xticklabel', {'Low','High'})

saveFigureIteration(f_persev, ...
                    '/Users/bib002/Documents/git-repositories/ed_perseveration/figures', ...
                    'theory_persev')