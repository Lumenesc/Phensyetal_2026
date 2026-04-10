%%%%%%%%%%%%%%%% plot_calcBustCorrelations.m %%%%%%%%%%%%%%%%%%%%
% Batch-Script to examine the correlation between LFP and GEVI signals
% during burst events.
%
% 

%% Initialize
tempoFreq = 40; 
range1 = 20;
lfpFreq = tempoFreq;
range2 = range1;

xCorrWin = 0.100; %Length of time (in seconds) to perform cross correlations on
maxLagtxl = 1/80; %Maximum lag (in seconds) for tempo x lfp correlation
stdDevThrsh = 2.0; %Threshold for how strong amplitude must be to count as a burst
minBurstDur = (1/40)*2; %Minimum required sustained burst activity (in seconds) to count as a burst

%% Calculate Burst Analysis for %%%%%%%%%%% PFC→MD Cohorts %%%%%%%%%%%
brstCorr = cell.empty(0,1); bursts = cell.empty(0,1);
% Get a list of all files and folders in this folder.
files = dir('C:\Users\AJPhe\Box\Postdoctoral Work\Gamma Oscillation Project\TEMPO Data\LFP x PFC-MD\Aligned TEMPO x LFP Data');
files(ismember( {files.name}, {'.', '..'})) = [];  %remove . and ..
for k = 1 : length(files) 
    
    load([files(k).folder '\' files(k).name]);
    fprintf('Performing burst analysis on data file #%d/%d = %s...\n', k, length(files), files(k).name);

    [brstCorr{k,1},~,bursts{k,1}] = calcBurstCorrelations(data,tempoFreq,range1,lfpFreq,range2,xCorrWin,maxLagtxl,stdDevThrsh,minBurstDur);
    
end

%% Calculate Burst Analysis for %%%%%%%%%%% PFC-PVI Cohorts %%%%%%%%%%%
brstCorr = cell.empty(0,1); brstCorr2 = cell.empty(0,1); bursts = cell.empty(0,1);
% Get a list of all files and folders in this folder.
files = dir('C:\Users\AJPhe\Box\Postdoctoral Work\Gamma Oscillation Project\TEMPO Data\LFP + PVI - Carl and Kathleen\carl_lfp+tempo\Day 1 Tanks');
files(ismember( {files.name}, {'.', '..'})) = [];
dirFlags = [files.isdir];
subFolders = files(dirFlags);

for k = 1 : length(files) 

    data = SEV2mat(strcat(subFolders(k).folder,'\',subFolders(k).name));
    data = combineLIAWav1(data);

    fprintf('Performing burst analysis on data file #%d/%d = %s...\n', k, length(files), files(k).name);

    [brstCorr{k,1},brstCorr2{k,1},bursts{k,1}] = calcBurstCorrelations(data,tempoFreq,range1,lfpFreq,range2,xCorrWin,maxLagtxl,stdDevThrsh,minBurstDur);
    
end


%% Combine brstCorrs and Plot Data

brstCorrAll = combinebrstCorr(brstCorr2);

%close all

figure
%corraxis=[-180:7.2:180];
corraxis=[-180:180/13:180];

for k=1:length(brstCorr)
subplot(length(brstCorr)+1,3,1+((k-1)*3))
bar(corraxis,mean(brstCorr{k,1}{2,1},2))
title('Gamma Bursts (30-50hz)')
ylabel(sprintf('Animal %d\nMean X-Corr',k))
xticks([-180:36:180])
ylim([-0.03 0.03])
subplot(length(brstCorr)+1,3,2+((k-1)*3))
bar(corraxis,mean(brstCorr{k,1}{2,2},2))
title('False Bursts (30-50hz)')
xticks([-180:36:180])
ylim([-0.03 0.03])
subplot(length(brstCorr)+1,3,3+((k-1)*3))
bar(corraxis,mean(brstCorr{k,1}{2,3},2))
title('Shuffled Data (30-50hz)')
xticks([-180:36:180])
ylim([-0.03 0.03])
end


subplot(length(brstCorr)+1,3,1+((k)*3))
bar(corraxis,mean(brstCorrAll{2,1},2))
title('Gamma Bursts (30-50hz)')
ylabel(sprintf('Average\nMean X-Corr'))
xticks([-180:36:180])
ylim([-0.02 0.02])
subplot(length(brstCorr)+1,3,2+((k)*3))
bar(corraxis,mean(brstCorrAll{2,2},2))
title('False Bursts (30-50hz)')
xticks([-180:36:180])
ylim([-0.02 0.02])
subplot(length(brstCorr)+1,3,3+((k)*3))
bar(corraxis,mean(brstCorrAll{2,3},2))
title('Shuffled Data (30-50hz)')
xticks([-180:36:180])
ylim([-0.02 0.02])
sgtitle(sprintf('LFP x TEMPO Gamma (30-50hz) Cross-Correlation\n %d Bursts: Min Duration: %dms. StdDev Thresh: %.1f',length(brstCorrAll{2,1}),round(minBurstDur*1000),stdDevThrsh))

%% Combine brstCorrs across different timepoints and Plot Mean Data for these different timepoints

%Apply behavioral timepoint constraints to burst arrays (e.g. only keep post-outcome RS bursts)
TSTART = 1; DIG = 2; OUTCOME = 3; TEND = 4; ITIEND = 5; CORRECT = 0; ERROR = 1;
trRange_labels = {'T.Start','DigStart','Outcome','T.End','ITI.END'};
CE_labels = {'Correct','Error'};

%plotset = [{'ALL',[TSTART ITIEND],[CORRECT ERROR]};{'IA',[TSTART DIG],[CORRECT ERROR]};{'IA',[OUTCOME TEND],[CORRECT]};{'IA',[OUTCOME TEND],[ERROR]};{'RS',[TSTART DIG],[CORRECT ERROR]};{'RS',[OUTCOME TEND],[CORRECT]};{'RS',[OUTCOME TEND],[ERROR]}];
%plotset = [{'IA',[OUTCOME TEND],[ERROR]};{'RS',[OUTCOME TEND],[ERROR]};{'IA',[TEND ITIEND],[ERROR]};{'RS',[TEND ITIEND],[ERROR]}];
plotset = [{'IA',[OUTCOME TEND],[ERROR]};{'RS',[OUTCOME TEND],[ERROR]}];

brstCorr_sel = brstCorr;

figure
corraxis=[-180:180/13:180];

burstCount = nan(size(brstCorr_sel,1),size(plotset,1));
burstLength = nan(size(brstCorr_sel,1),size(plotset,1));
brstCorr_Const = cell.empty(1,0);
brst_idx = cell.empty(1,0);
brstCorrAll = cell.empty(1,0);
for set = 1:size(plotset,1)
    taskphase = plotset{set,1}; %IA, RS, or ALL
    trRange = plotset{set,2};
    CEOPT =  plotset{set,3};
    [brstCorr_Const{1,set},trialCount,brst_idx{1,set}] = constrainbrstCorr(taskphase,trRange,CEOPT,timepoints,bursts,brstCorr_sel);
    brstCorrAll{1,set} = combinebrstCorr(brstCorr_Const{1,set});

    setTitle = sprintf('Phase: %s\n%s to %s\n%s\n(%d bursts)',taskphase,trRange_labels{trRange(1)},trRange_labels{trRange(2)},cell2mat(CE_labels(CEOPT+1)),size(brstCorrAll{1,set}{2,1},2));

    subplot(size(plotset,1),1,set)
    bar(corraxis,mean(brstCorrAll{1,set}{2,1},2))
    title('Gamma Bursts (30-50hz)')
    ylabel(sprintf('%s\nMean X-Corr',setTitle))
    xticks([-180:36:180])
    ylim([-0.05 0.05])

    %Grab the number of bursts per condition and mouse
    for k=1:size(brstCorr_Const{1,set},1)
        burstCount(k,set) = size(brstCorr_Const{1,set}{k,1}{2,1},2);
        burstLength(k,set) = mean(bursts{k,1}(brst_idx{1,set}{k,1},2)-bursts{k,1}(brst_idx{1,set}{k,1},1));
    end
end


sgtitle(sprintf('LFP x TEMPO Gamma (30-50hz) Cross-Correlation Across Different Behavioral Epochs\n Min Burst Duration: %dms. StdDev Thresh: %.1f',round(minBurstDur*1000),stdDevThrsh))




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUPPORT FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%% Function to combine TEMPO LIA and LFP Wav1 signals for carl/kathleen recorded data %%%%%%%
function data_out = combineLIAWav1(SEVdata)

    tempoData = SEVdata.x.data;
    lfpData = SEVdata.Wav1.data;

    if size(tempoData,2) > size(lfpData,2)
        fs_diff = round(size(tempoData,2)/size(lfpData,2));
        tempoData = downsample(tempoData',fs_diff)';
        fs = SEVdata.x.fs/fs_diff;
    else
        fs_diff = round(size(lfpData,2)/size(tempoData,2));
        lfpData = downsample(lfpData',fs_diff)';
        fs = SEVdata.Wav1.fs/fs_diff;
    end


    %Truncate data to same length
    if size(tempoData,2) > size(lfpData,2)
        tempoData = tempoData(:,1:size(lfpData,2));
    elseif size(lfpData,2) > size(tempoData,2)
        lfpData = lfpData(:,1:size(tempoData,2));
    end

    data_out = struct('lfp_data',lfpData,'tempo_data',tempoData,'fs',fs);

end


%% %%%%% Function to select bursts within constrained behavioral timepoints %%%%%%%
function [brstCorr_Const,trialCount,brst_idx] = constrainbrstCorr(taskphase,trRange,CEOPT,timepoints,bursts,brstCorr)
CLMN_IA = 2; CLMN_RS = 3; CLMN_CE = 6;
trRange_labels = {'T.Start','DigStart','Outcome','T.End','ITI.END'};
CE_labels = {'Correct','Error'};


tmpConsts = cell.empty(0,1);
%Build a set of constraints of timepoints
trialCount = zeros(size(timepoints,1));
for k = 1:size(timepoints,1)
    n=0;
    tmpConsts{k,1} = cell.empty(0,1);
    if strcmp(taskphase,'IA') || strcmp(taskphase,'ALL')
        for tr = 1:size(timepoints{k,CLMN_IA},1)
            if ismember(timepoints{k,CLMN_IA}{tr,CLMN_CE},CEOPT)
                n = n+1;
                tmpConsts{k,1}{n,1} = timepoints{k,CLMN_IA}{tr,trRange(1)}:timepoints{k,CLMN_IA}{tr,trRange(2)};
            end
        end
    end
    if strcmp(taskphase,'RS') || strcmp(taskphase,'ALL')
        for tr = 1:size(timepoints{k,CLMN_RS},1)
            if ismember(timepoints{k,CLMN_RS}{tr,CLMN_CE},CEOPT)
                n = n+1;
                tmpConsts{k,1}{n,1} = timepoints{k,CLMN_RS}{tr,trRange(1)}:timepoints{k,CLMN_RS}{tr,trRange(2)};
            end
        end
    end
    trialCount(k) = n;
end

brst_idx = cell.empty(0,1);
brstCount = [0 0];
for k=1:size(bursts,1)
    brst_idx{k,1} = false(size(bursts{k,1},1),1);
    for brst = 1:size(bursts{k,1},1)
        if any(cellfun(@(x) ismember(round(bursts{k,1}(brst,1)),x),tmpConsts{k,1}))
            brst_idx{k,1}(brst) = true;
        end
    end
    brstCount(1) = brstCount(1) + size(brst_idx{k,1},1);
    brstCount(2) = brstCount(2) + sum(brst_idx{k,1});
end

fprintf('%d out of %d bursts in taskphase: %s and trial range: %s through %s on %s trials\n',brstCount(2),brstCount(1),taskphase,trRange_labels{trRange(1)},trRange_labels{trRange(2)},cell2mat(CE_labels(CEOPT+1)));

brstCorr_Const = brstCorr;
for k = 1:size(brstCorr,1)
    for i=1:size(brstCorr{k,1},2)
        brstCorr_Const{k,1}{2,i} = brstCorr{k,1}{2,i}(:,brst_idx{k,1});
    end
end

end



%% %%%%% Function to combine all of the brstCorrs across mice %%%%%%%
function brstCorrAll = combinebrstCorr(brstCorr)

brstCorrAll = cell(2,3);
brstCorrAll{1,1} = 'Bursts'; brstCorrAll{1,2} = 'False Bursts'; brstCorrAll{1,3} = 'Shuffled';
for k= 1 :length(brstCorr)
    if k == 1
        for i=1:3
            brstCorrAll{2,i} = brstCorr{k,1}{2,i};
        end
    else
        for i=1:3
            brstCorrAll{2,i} = [brstCorrAll{2,i} brstCorr{k,1}{2,i}];
        end
    end
end

end
%EOF