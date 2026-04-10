%{

AnalyzeTEMPO_RS_Data v1

This script consolidates the SigMatrix data set acquired through
plot_temposcript_XCORR and the timepoints data set acquired through
loadTimepoints. Please have these two files ready prior to running script.

This version utilizes PHASE SHIFTED SigMatrix data - e.g.
temposcript_XCORR_PS was used to generate the SigMatrix

This version only analyzes TEMPO with a single voltage sensor in two sites, 
and is not compatible with 2C-TEMPO data sets. If you would like to analyze
each voltage signal from a  2C-TEMPO dataset separately, please first use 
split2C.m to convert the 2C-SigMatrix to 2 separate 1C-SigMatrices.

This version is compatible with - but does not require - the trial
classifications designed by Caitriona Costello - however the program will
expect you to have these datasets already prepared - if you do not have
them simply decline to load.

Description of Analysis: All SigMatrix values consist of correlation value
between recording sites calculated across a number of phase shifts. This
code first reduces this set to the 'peak' phase + r-value for each
timepoint.  It then restricts analysis to a user-provided phase range (e.g.
-45 : 45 degrees) therefore quantifying the phase-specific synchrony.
Finally the code averages this synchrony within specific time-frames
defined by the RS task structure and the timepoints data set to allow the
user to compare synchrony between task-relevant epochs. 

Modifiable parameters: There are numerous parameters that can be changed by
the user, therefore after each run there will be a parameters file automatically
exported with a timestamp to allow the user to refer back to that specific
run - in addition, all graphs generated have a timestamp to allow for
cross-referencing to the parameters file

Data Interaction: After populating the initial graphs, the user can play
with the data through a set of auxilary scripts, as I complete these
scripts I will provide a brief description here:

---

---


Aarron Phensy - Sohal Lab - November 2024




%}

%% Initialize
for i=1:15
    for j=1:i
        fprintf('.');
    end
    fprintf('\n');
end
fprintf('Executing AnalyzeTEMPO_RS_Functions\n\n');
fprintf(['This script will help you build both behavioral IDX and synchrony datasets needed to investigate synchrony at\n' ...
    'specific behavioral epochs - for example investigate synchrony in correct versus error outcomes during Rule Shift. After\n' ...
    'generating these datasets you will be able to plot graphs using plotSynchxIDX.m.\n'])
fprintf('\nPress any key to get started.\n');
pause();
clear i j

addpath('AnalyzeTEMPO_RS_Functions')

%% Load user TEMPO and RS data

if exist('cohortData','var')
    pendingAnswer = true;
    while pendingAnswer
        userInput = input('Load new user data? (Y/N): ',"s");
        if strcmp(userInput,'N') || strcmp(userInput,'n')
            fprintf('Continuing with existing user data...\n\n');
            pendingAnswer = false;
        elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
            cohortData = loadUserTEMPORSData();
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
else
    cohortData = loadUserTEMPORSData();
end

if exist('trialSectionTimes','var')
    pendingAnswer = true;
    while pendingAnswer
        userInput = input('Redefine trial section ranges? (Y/N): ',"s");
        if strcmp(userInput,'N') || strcmp(userInput,'n')
            fprintf('Continuing with existing trial section ranges...\n\n');
            pendingAnswer = false;
        elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
            trialSectionTimes = setTrialSections();
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
else
    trialSectionTimes = setTrialSections();
end

clear userInput pendingAnswer

%% Configurations for TEMPO Phase-Locked Synchrony Analysis
fprintf('\n\nNow to perform the synchrony analysis time-locked to a given phase range.\n');


if exist('phaseset','var')
    pendingAnswer = true;
    while pendingAnswer
        userInput = input('Redefine phase range(s)? (Y/N): ',"s");
        if strcmp(userInput,'N') || strcmp(userInput,'n')
            fprintf('Continuing with existing phase range(s)...\n\n');
            pendingAnswer = false;
        elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
            numPhases = input('How many phase ranges do you want to analyze? E.g. [-45,45] is a single range for 0-Lag analysis.  ');
            phaseset = cell(1,0);
            for i=1:numPhases
                fprintf('\nSetting ranges for phase set # %d\n',i);
                phaseset{i}(1) = input('Please enter the lower bound of the phase range (in degrees): ');
                phaseset{i}(2) = input('Please enter the upper bound of the phase range (in degrees): ');
            end
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
else
    numPhases = input('How many phase ranges do you want to analyze? E.g. [-45,45] is a single range for 0-Lag analysis.  ');
    phaseset = cell(1,0);
    for i=1:numPhases
        fprintf('\nSetting ranges for phase set # %d\n',i);
        phaseset{i}(1) = input('Please enter the lower bound of the phase range (in degrees): ');
        phaseset{i}(2) = input('Please enter the upper bound of the phase range (in degrees): ');
    end
end
fprintf('Synchrony will be performed on the following phase ranges: ')
for i=1:length(phaseset)
    fprintf('[%d, %d] ',phaseset{i}(1),phaseset{i}(2))
end
fprintf('\nPress any key to continue.\n')
pause();

pendingAnswer = true;
while pendingAnswer
    userInput = input('\n\nNormalize data to baseline? (Y/N): ',"s");
    if strcmp(userInput,'N') || strcmp(userInput,'n')
        fprintf('Normalization Off - Averaging r-values.\n');
        normON = false;
        pendingAnswer = false;
    elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
        fprintf('Normalizating data to baseline.\n');
        normON = true;
        pendingAnswer = false;
    else
        fprintf('Invalid Input\n');
    end
end

pendingAnswer = true;
while pendingAnswer
    userInput = input('\n\nCalculate Synchrony to Reference-Channel? Warning, doubles the # of plots. (Y/N): ',"s");
    if strcmp(userInput,'N') || strcmp(userInput,'n')
        fprintf('Reference-Channel Synchrony Off.\n');
        calculateTDTomato = false;
        pendingAnswer = false;
    elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
        fprintf('Signal to Reference-Channel Synchrony will be added to cohorts.\n');
        calculateTDTomato = true;
        pendingAnswer = false;
    else
        fprintf('Invalid Input\n');
    end
end

% Peform TEMPO Phase-Locked Synchrony Analysis
analysisCounter=0;
synchronyData = cell(2,1);
behavioralIDXData = cell(2,1);

%metadata_SIGSHF = cell(size(cohortData),1);
for cohort=1:size(cohortData,1)

    [IDX_set,metadata_IDX] = getBehavioralIndices(cohortData{cohort,3});
    behavioralIDXData{1,cohort} = sprintf('%s IDX_Set',cohortData{cohort,1});
    behavioralIDXData{2,cohort} = IDX_set;
    behavioralIDXData{3,cohort} = {'PrevNonConfPost','PrevConfPost','Pre','Dig','Out','Post','ITISt','ITIEnd','NxtPre','NxtNCPre','NxtCPre','Trial'}; %In a future release, make this customizable..

    for i=1:length(phaseset)
        phaserange = phaseset{i};
        
        analysisCounter = analysisCounter+1;
        
        [synchronyValues] = calcTEMPOSynch_PhaseLocked(cohortData{cohort,4},cohortData{cohort,3},phaserange,trialSectionTimes,normON);
        synchronyData{1,cohort} = sprintf('%s',cohortData{cohort,1});
        synchronyData{2,cohort}{1,i} = sprintf('%d:%d',phaserange(1),phaserange(2));
        synchronyData{2,cohort}(2:length(synchronyValues)+1,i) = synchronyValues';
    
    end

    if calculateTDTomato == true
        behavioralIDXData{1,cohort+size(cohortData,1)} = sprintf('%s IDX_Set',cohortData{cohort,1});
        behavioralIDXData{2,cohort+size(cohortData,1)} = IDX_set;
        behavioralIDXData{3,cohort+size(cohortData,1)} = {'PrevNonConfPost','PrevConfPost','Pre','Dig','Out','Post','ITISt','ITIEnd','NxtPre','NxtNCPre','NxtCPre','Trial'}; %In a future release, make this customizable..
        for i=1:length(phaseset)
            phaserange = phaseset{i};

            analysisCounter = analysisCounter+1;

            referenceData = cohortData{cohort,4};
            referenceData(:,[2,3]) = referenceData(:,[3,2]); %Reference data just flips the columns used for analysis

            [synchronyValues] = calcTEMPOSynch_PhaseLocked(referenceData,cohortData{cohort,3},phaserange,trialSectionTimes,normON);
            synchronyData{1,cohort+size(cohortData,1)} = sprintf('%s-REF',cohortData{cohort,1});
            synchronyData{2,cohort+size(cohortData,1)}{1,i} = sprintf('%d:%d',phaserange(1),phaserange(2));
            synchronyData{2,cohort+size(cohortData,1)}(2:length(synchronyValues)+1,i) = synchronyValues';
        end

    end
    
    % metadata_SIGSHF{j,1} = cohortData{cohort,1};
    % metadata_SIGSHF{j,2} = sigSHFCompComplex{1,4};

end

% metadata_cohortInfo = struct('cohortName',cohortName,'cohortset',{cohortset},'phaseset',{phaseset},'trsectiontimes',{trSecTimes},'NormalizeSynchValues',normON);
% metadata_All = struct('Cohort_Meta',metadata_cohortInfo,'SigSHF_Meta',metadata_SIGSHF,'IDX_Meta',metadata_IDX);
% 
% if strcmp(saveMeta,'on')
%     t = datetime;
%     t.Format = "yyyy-MM-dd-hh-mm-ss";
%     save(sprintf('metadata\\%s_%s_calcTaskSingle',t,cohortName),'metadata_All');
% end
% 
% OutSynchFinalSet = {synchronyData,OutSynchFinalCombVectors,{},IDX_set_cohort,OUT_behaviordata_cohort,metadata_All};
% %assignin('base',sprintf('OutSynchFinalSet%s',cohortName),OutSynchFinalSet);
% calcTaskPhaseSigFunctionsCombCohorts

clear metadata_All metadata_cohortInfo metadata_SIGSHF metadata_IDX 
clear i cohort IDX_set synchronyValues normON analysisCounter calculateTDTomato phaserange numPhases
clear userInput pendingAnswer referenceData