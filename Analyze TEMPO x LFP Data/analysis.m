%% Get Timepoints Constraints to apply to sigMatrix
sigMatrix_Z = sigMatrix_Z_bkp;
sigSHFComp = sigSHFComp_bkp;

trRange_labels = {'T.Start','DigStart','Outcome','T.End','ITI.END'};
CE_labels = {'Correct','Error'};

taskphase = 'RS'; %IA, RS, or ALL
trRange = [3 4];
CEOPT = [0];

tmpConsts = getTPConstraints(taskphase,trRange,CEOPT,timepoints);

for k=1:size(tmpConsts,1)
    tmpConsts{k,1} = round(tmpConsts{k,1}*(window/102));
end

for k=1:size(sigMatrix,1)
    for i=3:size(sigMatrix,2)
        sigMatrix_Z{k,i} = sigMatrix_Z{k,i}(:,tmpConsts{k,1},:);
    end
    for j=1:size(sigSHFComp,2)
        sigSHFComp{k,j} = sigSHFComp{k,j}(:,tmpConsts{k,1});
    end
end

%%
TESTARRAY = LFPAECARRAY;
t=t+1;
p(t) = ranksum(TESTARRAY(S), TESTARRAY(~S));
means(t,1) = mean(TESTARRAY(S));
means(t,2) = mean(TESTARRAY(~S));

obs = median(TESTARRAY(S)) - median(TESTARRAY(~S));
nPerm = 10000;
null = zeros(nPerm,1);
for n = 1:nPerm
    Sper = S(randperm(numel(S)));
    null(n) = median(TESTARRAY(Sper)) - median(TESTARRAY(~Sper));
end
p_perm(t) = mean(abs(null) >= abs(obs));

fprintf('LFP Phase-Locking Value: Mean PLV During Synch: %f vs ~Synch: %f, Ranksum p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_perm(t));
figure; hold on

% Histogram of empirical null distribution
histogram(null, 40, ...
    'Normalization','probability', ...
    'FaceColor',[0.7 0.7 0.7], ...
    'EdgeColor','none');

% Vertical line for observed statistic
yl = ylim;
plot([obs obs], yl, 'r-', 'LineWidth', 2);

xlabel('\Delta PLV (median Sync - median NonSync)')
ylabel('Probability')
title(sprintf('Empirical Null Distribution (perm p = %.4f)', p_perm(t)))

box off
set(gca,'FontSize',12)

%% Perform Ranksum and Label Permutation Tests across Measures
fprintf('\nPerforming Ranksum and Label Permutation Test across measures \n');
clear p p_perm means
S = logical.empty(1,0);
POWERARRAY = double.empty(1,0);
LFPxTEMPOARRAY = double.empty(1,0);
LFPxCorrARRAY = double.empty(1,0);
LFPPLVARRAY = double.empty(1,0);
LFPAECARRAY = double.empty(1,0);
LFPwPLIARRAY = double.empty(1,0);
LFPxCorrRawARRAY = double.empty(1,0);
LFPxCorrARRAYLagsARRAY = double.empty(1,0);
for k = animalSet
    S = [S sigSHFComp{k,1}==1;];
    POWERARRAY = [POWERARRAY mean(sigMatrix_Z{k,6},1)];
    LFPxTEMPOARRAY = [LFPxTEMPOARRAY mean(sigMatrix_Z{k,5},1)];
    LFPxCorrARRAY = [LFPxCorrARRAY max(sigMatrix_Z{k,4}(:,:,41),[],3)];
    LFPPLVARRAY = [LFPPLVARRAY sigMatrix_Z{k,8}];
    LFPAECARRAY = [LFPAECARRAY sigMatrix_Z{k,9}];
    LFPwPLIARRAY = [LFPwPLIARRAY sigMatrix_Z{k,12}];
    LFPxCorrRawARRAY = [LFPxCorrRawARRAY max(sigMatrix_Z{k,11}(:,:,41),[],3)];
end

t=0;
t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(POWERARRAY,S);
fprintf('LFP Power: Mean Power During Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong> , KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPxTEMPOARRAY,S);
fprintf('LFPXTEMPO XCORR: Mean Correlation During Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPxCorrARRAY,S);
fprintf('LFPXLFP XCORR: Mean Correlation During Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPPLVARRAY,S);
fprintf('LFP Phase-Locking Value: Mean PLV During Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPAECARRAY,S);
fprintf('LFP AEC: Mean AEC During Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPwPLIARRAY,S);
fprintf('LFP wPLI: Mean wPLI During Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPxCorrRawARRAY,S);
fprintf('LFPxLFP xCorr (Full-Band) Value: Mean Correlation During Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

% figure; hold on
% histogram(lag_ms(~S), -25:1:25, 'Normalization','probability');
% histogram(lag_ms(S),  -25:1:25, 'Normalization','probability');
% xlabel('Lag (ms)'); ylabel('Probability');
% legend('~Sync','Sync');
% title('Lag at peak LFP-LFP correlation');

%  Perform Ranksum and Label Permutation Tests for TEMPO X-Corr when significantly high LFP X-Corr
fprintf('\nPerforming Ranksum and Label Permutation Test across measures when LFP has high gamma-synch\n');
clear p p_perm
S = logical.empty(1,0);
TEMPOxCorrARRAY = double.empty(1,0);
LFPxTEMPOARRAY2 = double.empty(1,0);
for k = animalSet
    S = [S sigSHFComp{k,3}==1;];
    TEMPOxCorrARRAY = [TEMPOxCorrARRAY sigMatrix_Z{k,3}(1,:)];
    LFPxTEMPOARRAY2 = [LFPxTEMPOARRAY2 mean(sigMatrix_Z{k,5},1)];
end

t=0;

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(TEMPOxCorrARRAY,S);
fprintf('TEMPO Cross-Correlation: Mean Correlation During LFP Gamma-Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, Perm p = %f \n',means(t,1),means(t,2),p(t),p_perm(t));

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPxTEMPOARRAY2,S);
fprintf('LFPXTEMPO XCORR: Mean Correlation During LFP-Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

%  Perform Ranksum and Label Permutation Tests for TEMPO X-Corr when significantly high LFP AEC - e.g. when gamma bursting is correlated in L/R LFPs
fprintf('\nPerforming Ranksum and Label Permutation Test across measures - High AEC \n');
clear p p_perm
S = logical.empty(1,0);
TEMPOxCorrARRAYAEC = double.empty(1,0);
LFPxTEMPOARRAYAEC = double.empty(1,0);
for k = animalSet
    S = [S sigSHFComp{k,4}==1;];
    TEMPOxCorrARRAYAEC = [TEMPOxCorrARRAYAEC sigMatrix_Z{k,3}(1,:)];
    LFPxTEMPOARRAYAEC = [LFPxTEMPOARRAYAEC mean(sigMatrix_Z{k,5},1)];
end

t=0;

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(TEMPOxCorrARRAYAEC,S);
fprintf('TEMPO Cross-Correlation: Mean Correlation During LFP AEC-Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, Perm p = %f \n',means(t,1),means(t,2),p(t),p_perm(t));

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPxTEMPOARRAYAEC,S);
fprintf('LFPXTEMPO XCORR: Mean Correlation During AEC-Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

%  Perform Ranksum and Label Permutation Tests for TEMPO X-Corr when gamma bursts exist in (concatenating L and R)
fprintf('\nPerforming Ranksum and Label Permutation Test across measures (Bursts L/R concatenated)\n');
clear p p_perm
S = logical.empty(1,0);
TEMPOxCorrARRAYBURST = double.empty(1,0);
LFPxTEMPOARRAYBURST = double.empty(1,0);
for k = animalSet
    S = [S sigSHFComp{k,5}(1,:) sigSHFComp{k,5}(2,:)];
    LFPxTEMPOARRAYBURST = [LFPxTEMPOARRAYBURST sigMatrix_Z{k,5}(1,:) sigMatrix_Z{k,5}(2,:)];
end

t=0;


t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPxTEMPOARRAYBURST,S);
fprintf('LFPXTEMPO XCORR: Mean Correlation During LFP Gamma Burst: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));

%  Perform Ranksum and Label Permutation Tests for TEMPO X-Corr when gamma bursts exist in (combining L and R via OR)
fprintf('\nPerforming Ranksum and Label Permutation Test across measures (Bursts L/R combine via OR operation) \n');
clear p p_perm
S = logical.empty(1,0);
TEMPOxCorrARRAYBURST2 = double.empty(1,0);
LFPxTEMPOARRAYBURST2 = double.empty(1,0);
for k = animalSet
    S = [S or(sigSHFComp{k,5}(1,:),sigSHFComp{k,5}(2,:))];
    TEMPOxCorrARRAYBURST2 = [TEMPOxCorrARRAYBURST2 sigMatrix_Z{k,3}(1,:)];
    LFPxTEMPOARRAYBURST2 = [LFPxTEMPOARRAYBURST2 mean(sigMatrix_Z{k,5},1)];
end

t=0;

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(TEMPOxCorrARRAYBURST2,S);
fprintf('TEMPO Cross-Correlation: Mean Correlation During LFP AEC-Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, Perm p = %f \n',means(t,1),means(t,2),p(t),p_perm(t));

t=t+1;
[p(t),p_KS(t),p_perm(t),means(t,:)] = perform_stats(LFPxTEMPOARRAYBURST2,S);
fprintf('LFPXTEMPO XCORR: Mean Correlation During AEC-Synch: %f vs ~Synch: %f, <strong> Ranksum p = %f</strong>, KS p = %f, Perm p = %f \n',means(t,1),means(t,2),p(t),p_KS(t),p_perm(t));


%% Support Functions

function [p,p_KS,p_perm,means] = perform_stats(data,S)

TESTARRAY = data;
p = ranksum(TESTARRAY(S), TESTARRAY(~S));
means(1,1) = mean(TESTARRAY(S));
means(1,2) = mean(TESTARRAY(~S));

obs = median(TESTARRAY(S)) - median(TESTARRAY(~S));
nPerm = 10000;
null = zeros(nPerm,1);
for n = 1:nPerm
    Sper = S(randperm(numel(S)));
    null(n) = median(TESTARRAY(Sper)) - median(TESTARRAY(~Sper));
end
p_perm = mean(abs(null) >= abs(obs));

[~,p_KS,~] = kstest2(TESTARRAY(S), TESTARRAY(~S));

end


function tmpConsts = getTPConstraints(taskphase,trRange,CEOPT,timepoints)

CLMN_IA = 2; CLMN_RS = 3; CLMN_CE = 6;
trRange_labels = {'T.Start','DigStart','Outcome','T.End','ITI.END'};
CE_labels = {'Correct','Error'};


tmpConsts = cell.empty(0,1);
%Build a set of constraints of timepoints
n=0;
for k = 1:size(timepoints,1)
    tmpConsts{k,1} = double.empty(0,1);
    if strcmp(taskphase,'IA') || strcmp(taskphase,'ALL')
        for tr = 1:size(timepoints{k,CLMN_IA},1)
            if ismember(timepoints{k,CLMN_IA}{tr,CLMN_CE},CEOPT)
                n = n+1;
                tmpConsts{k,1} = [tmpConsts{k,1},timepoints{k,CLMN_IA}{tr,trRange(1)}:timepoints{k,CLMN_IA}{tr,trRange(2)}];
            end
        end
    end
    if strcmp(taskphase,'RS') || strcmp(taskphase,'ALL')
        for tr = 1:size(timepoints{k,CLMN_RS},1)
            if ismember(timepoints{k,CLMN_RS}{tr,CLMN_CE},CEOPT)
                n = n+1;
                tmpConsts{k,1} = [tmpConsts{k,1},timepoints{k,CLMN_RS}{tr,trRange(1)}:timepoints{k,CLMN_RS}{tr,trRange(2)}];
            end
        end
    end
end

end