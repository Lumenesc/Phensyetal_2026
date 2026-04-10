%%
%STEP 1
%Initialize when first opening matlab
addpath('distFig v5.5')
addpath('pierremegevand-watsons_u2-cce8e55')
%%

%STEP 2
%First - select folder to populate each animal's session length (in
%samples)
files = dir('C:\Users\AJPhe\Box\Postdoctoral Work\Gamma Oscillation Project\TEMPO Data\TCM TEMPO Files\RunMe');
files(ismember( {files.name}, {'.', '..'})) = [];  %remove . and ..
dirFlags = [files.isdir]; % Get a logical vector that tells which is a directory.
subFolders = files(dirFlags); % Extract only those that are directories.

%Collect session length for each subdirectory folder ##NOTE: All subdirectories must be SEV container folders
clear sessionlength
for k = 1 : length(subFolders)
    [data] = SEV2mat(strcat(subFolders(k).folder,'\',subFolders(k).name));
    fprintf('Collecting data on sub folder #%d/%d = %s...\n', k, length(subFolders), subFolders(k).name);
    if k==1
        sessionlength(k,1) = length(data.x.data(1,:))+1; %%%%%%%%%NOTE: Need to add 1 sample to the length of the first animal due to how 'masterx' was built - needs to be corrected
    else
        sessionlength(k,1) = length(data.x.data(1,:));
    end
end
[data] = SEV2mat(strcat(subFolders(1).folder,'\',subFolders(1).name));
clear files dirFlags subFolders k

%%

%STEP 3

MF_fs = data.x.fs;
MF_inc = win; %in seconds
MF_win = inc;
MF_dataOffset = ceil(MF_win/MF_inc);
MF_incsize = round(MF_inc*MF_fs);
MF_npts = sum(sessionlength(:,1)); 

%Load timepoints - padITIEnd and convert to inc (e.g. 250ms) w/ convTimepointFS(timepoints,fs)
[file,path] = uigetfile({'*.mat'},'Timepoints Struct'); 
load([path file]);
timepoints = padITIEnd(timepoints,10); %Pad 10 seconds to missing ITI ends - note, if recordings are terminated too close to final trial end, you may need to shorten this.
timepoints = convTimepointFS(timepoints,1/MF_inc);

%Load IDX_set (populated from calcTaskPhaseSigMinimal)
[file,path] = uigetfile({'*.mat'},'IDX_set'); 
load([path file]);

for k=1:size(sessionlength,1)
    if k<size(sessionlength,1)
        sessionlength(k,2)=floor(sum(sessionlength(1:k),2)/MF_incsize);
    else
        sessionlength(k,2)=floor(sum(sessionlength(1:k),2)/MF_incsize)-MF_dataOffset;
    end

    if k==1
        sessionlength(k,3)=sessionlength(k,2);
    else
        sessionlength(k,3) = sessionlength(k,2)-sessionlength(k-1,2);
    end
end
   
clear file path data k 
%%

%Step 4

%IDX for trial-based info

%Trial selection (manual) - each [] is the set of trials to include
%Note - TCM3-2 is excluded in the main figure plots

clear trSecTimes
TSTART = 1; DIG = 2; OUTCOME = 3; TEND = 4; ITIEND = 5;
trSecTimes(1,:) = {'Pre-Decision','Dig','Outcome','Post-Outcome','ITI Start','ITI End','Max-PreDec'};
trSecTimes(2,:) = {[TSTART,0,DIG,-1],[DIG,-1,OUTCOME,-1],[OUTCOME,-2,OUTCOME,2],[OUTCOME,-1,TEND,0],[TEND,0,TEND,15],[ITIEND,-15,ITIEND,0],10}; %Time ranges (in seconds) for: Pre-Decision, Dig, Outcome, Post-Outcome, ITI-Start, ITI-End

%SET RANGES E.G. 1=Trial Start, 2=Dig, 3=Outcome, 4=Trial End, 5=ITI End.
tRange1 = 1;
tRange2 = 5;

divSplitLimit = [99 99]; %Maximum trials [before after] to include in the Divergent Split metric
divSplitOmit = [0 0]; %Trials to omit from the divergent trial [before after], incorporated to not let the trials right around the div contaminate comparisons of pre and post

IASplitLimit = [4 4]; %Maximum trials [start end] of the IA-Policy to include in the 'IA-Policy' split - note: this is during RS only
IASplitMinDist = 4; %Minimum allowable distance from start and end trials, Note: this may exclude animals who leave IA in less trials than the MinDist

clear tSelectionNames tSelection

%Final -3:IA:+5
tSelectionNames = {'Non-Conflict','Conflict'}; %Conflict vs Non-Conflict
tSelection{1} = {[2,3,5,8],[16,17,19,21,24],[9,12,14,16],[12,13,16],[4,6,8,10,11],[26,27,30,31],[2,3,6,8],[],[]};  
tSelection{2} = {[1,4,6,7],[18,20,22,23],[10,11,13,15],[11,14,15,17],[5,7,9,12],[25,28,29,32],[1,4,5,7],[],[]};  


clear IDX_time_subject %Subject ID (k)
clear IDX_time_task %BL == 1, IA == 2, RS == 3, Post-RS == 4
clear IDX_time_trial %Both IA and RS trial #s T.Start : ITI.End
clear IDX_time_phase %Non-TaskPhase == 0, %Pre-Dec == 1, Post-Out == 2, ITI == 3
clear IDX_time_CE %Non-TaskPhase == 0, Correct == 1, Error == 2
clear IDX_time_pCE %Non-TaskPhase == 0, p.Correct == 1, p.Error == 2
clear IDX_time_pair %Non-TaskPhase == 0, ReinfPair == 1, NewPair == 2
clear IDX_time_ppair %Non-TaskPhase/NoPrevTrial == 0, ReinfPair == 1, NewPair == 2
clear IDX_time_caitComb
clear IDX_time_caitTrans %Caitriona Classification Transitions - IA-Policy = 1; IA-Divergence = 2; Exploratory = 3; RS-Acquisition = 4; RS-Policy = 5
clear IDX_time_caitCombxPair
clear IDX_time_selectTrials
clear IDX_time_divSplit %Split trials as Pre-Divergent == 1, Divergent == 2, and Post-Divergent == 3
clear IDX_time_IASplit %Split trials as IA-Policy Start == 1, IA-Policy End == 2, and Divergent == 3
clear IDX_time_All

for k=1:size(sessionlength,1)

    %%IDX_time_All - All Timepoints
    IDX_time_All{k,1} = ones(1,sessionlength(k,3));

    %%IDX_time_selectTrials - Manual Selections
    IDX_time_selectTrials{k,1} = zeros(1,sessionlength(k,3));
    % for i=1:size(tSelection,2)
    %     for tr=tSelection{i}{k} %RS Only for now...
    %         timestamps = cell2mat(timepoints{k,3}(tr,1:5)); %Trial Start:ITI END
    %         IDX_time_selectTrials{k,1}(timestamps(1):timestamps(5)) = i;
    %     end
    % end

    %%IDX_time_subject - Subject ID (k)
    IDX_time_subject{k,1} = ones(1,sessionlength(k,3))*k;

    %%IDX_time_task - BL, IA, RS, Post-RS
    clear timestamps
    timestamps(1) = cell2mat(timepoints{k,2}(1,1)); %IAStart
    timestamps(2) = cell2mat(timepoints{k,3}(1,1)); %RSStart
    timestamps(3) = cell2mat(timepoints{k,3}(end,5)); %RSEnd
    IDX_time_task{k,1} = ones(1,sessionlength(k,3));
    IDX_time_task{k,1}(timestamps(1):timestamps(2)) = 2;
    IDX_time_task{k,1}(timestamps(2):timestamps(3)) = 3;
    IDX_time_task{k,1}(timestamps(3):end) = 4;

    %%IDX_time_trial - Trial # 1:nTrials, T.Start to ITI.End
    clear timestamps
    IDX_time_trial{k,1} = zeros(1,sessionlength(k,3));
    for tr=1:size(timepoints{k,2},1) %IA First
        for i=1:5
            timestamps(i) = cell2mat(timepoints{k,2}(tr,i));
        end
        IDX_time_trial{k,1}(timestamps(1):timestamps(5)) = tr;
    end
    for tr=1:size(timepoints{k,3},1) %RS Second
        for i=1:5
            timestamps(i) = cell2mat(timepoints{k,3}(tr,i));
        end
        IDX_time_trial{k,1}(timestamps(1):timestamps(5)) = tr+size(timepoints{k,2},1); %Trials keep counting from IA
    end
    
    %%IDX_time_phase - Non-TaskPhase == 0, %Pre-Dec == 1, Post-Out == 2, ITI == 3
    clear timestamps
    IDX_time_phase{k,1} = zeros(1,sessionlength(k,3));
    %IDX_time_phase{k,2} = zeros(1,sessionlength(k,3));
    for tr=1:size(timepoints{k,2},1) %IA First
        
        timestamps(:) = cell2mat(timepoints{k,2}(tr,1:5));

        if timestamps(2)-timestamps(1) > trSecTimes{2,7}/MF_inc %If Tstart:Dig is greater than max pre-decision time
            IDX_time_phase{k,1}(timestamps(2)-trSecTimes{2,7}/MF_inc:timestamps(trSecTimes{2,1}(3))+(trSecTimes{2,1}(4)/MF_inc)) = 1; %Just take the max pre-decision time from dig
        else
            IDX_time_phase{k,1}(timestamps(trSecTimes{2,1}(1))+(trSecTimes{2,1}(2)/MF_inc):timestamps(trSecTimes{2,1}(3))+(trSecTimes{2,1}(4)/MF_inc)) = 1; %Pre-Decision from t.start
        end
        IDX_time_phase{k,1}(timestamps(trSecTimes{2,4}(1))+(trSecTimes{2,4}(2)/MF_inc):timestamps(trSecTimes{2,4}(3))+(trSecTimes{2,4}(4)/MF_inc)) = 2; %Post-Outcome
        IDX_time_phase{k,1}(timestamps(4)+3:timestamps(5)) = 3;

    end
    for tr=1:size(timepoints{k,3},1) %RS Second
        
        timestamps(:) = cell2mat(timepoints{k,3}(tr,1:5));
        
        if timestamps(2)-timestamps(1) > trSecTimes{2,7}/MF_inc %If Tstart:Dig is greater than max pre-decision time
            IDX_time_phase{k,1}(timestamps(2)-trSecTimes{2,7}/MF_inc:timestamps(trSecTimes{2,1}(3))+(trSecTimes{2,1}(4)/MF_inc)) = 1; %Just take the max pre-decision time from dig
        else
            IDX_time_phase{k,1}(timestamps(trSecTimes{2,1}(1))+(trSecTimes{2,1}(2)/MF_inc):timestamps(trSecTimes{2,1}(3))+(trSecTimes{2,1}(4)/MF_inc)) = 1; %Pre-Decision from t.start
        end
        IDX_time_phase{k,1}(timestamps(trSecTimes{2,4}(1))+(trSecTimes{2,4}(2)/MF_inc):timestamps(trSecTimes{2,4}(3))+(trSecTimes{2,4}(4)/MF_inc)) = 2; %Post-Outcome
        IDX_time_phase{k,1}(timestamps(4)+3:timestamps(5)) = 3; %ITI

    end
    
    clear timestamps
    IDX_time_CE{k,1} = zeros(1,sessionlength(k,3));
    %IDX_time_CE{k,2} = zeros(1,sessionlength(k,3));
    IDX_time_pair{k,1} = zeros(1,sessionlength(k,3));
    IDX_time_caitComb{k,1} = zeros(1,sessionlength(k,3));
    IDX_time_caitTrans{k,1} = zeros(1,sessionlength(k,3));
    %IDX_time_caitComb{k,2} = zeros(1,sessionlength(k,3));
    IDX_time_caitCombxPair{k,1} = zeros(1,sessionlength(k,3));
    IDX_time_divSplit{k,1} = zeros(1,sessionlength(k,3));
    IDX_time_IASplit{k,1} =  zeros(1,sessionlength(k,3));
    for tr=1:size(timepoints{k,2},1) %IA
        timestamps = cell2mat(timepoints{k,2}(tr,1:5)); %Trial Start:ITI END
        if timepoints{k,2}{tr,6} == 0
            IDX_time_CE{k,1}(timestamps(1):timestamps(5)) = 1;
        elseif timepoints{k,2}{tr,6} == 1
            IDX_time_CE{k,1}(timestamps(1):timestamps(5)) = 2;
        end

        %IDX_time_caitComb - Non-Phase = 0; IA-Policy = 1; Divergent(S+E) =2; Exploratory = 3; RS-Policy = 4;
        IDX_time_caitComb{k,1}(timestamps(1):timestamps(5)) = IDX_set{2,4}{k,1}(tr);
        
        %Caitriona Classification Transitions - IA-Policy = 1; IA-Divergence = 2; Exploratory = 3; RS-Acquisition = 4; RS-Policy = 5
        if ~isnan(IDX_set{2,5}{k,1}(tr))
            IDX_time_caitTrans{k,1}(timestamps(1):timestamps(5)) = IDX_set{2,5}{k,1}(tr);
        end

    end
    for tr=1:size(timepoints{k,3},1) %RS
        timestamps = cell2mat(timepoints{k,3}(tr,1:5)); %Trial Start:ITI END
        if timepoints{k,3}{tr,6} == 0
            IDX_time_CE{k,1}(timestamps(tRange1):timestamps(tRange2)) = 1;
        elseif timepoints{k,3}{tr,6} == 1
            IDX_time_CE{k,1}(timestamps(tRange1):timestamps(tRange2)) = 2;
        end
        
        %IDX_time_pair - Non-Phase = 0; Non-Conflict = 1; Conflict = 2; 
        IDX_time_pair{k,1}(timestamps(tRange1):timestamps(tRange2)) = IDX_set{2,16}{k,1}(tr);

        %IDX_time_caitComb - Non-Phase = 0; IA-Policy = 1; Divergent(S+E) =2; Exploratory = 3; RS-Policy = 4;
        IDX_time_caitComb{k,1}(timestamps(tRange1):timestamps(tRange2)) = IDX_set{2,4}{k,2}(tr);

        %Caitriona Classification Transitions - IA-Policy = 1; IA-Divergence = 2; Exploratory = 3; RS-Acquisition = 4; RS-Policy = 5
        if ~isnan(IDX_set{2,5}{k,2}(tr))
            IDX_time_caitTrans{k,1}(timestamps(tRange1):timestamps(tRange2)) = IDX_set{2,5}{k,2}(tr);
        end

        %Divergent Split Trials - Using CaitTransition
        divTr = find(IDX_set{2,5}{k,2}==2,1,'first');
        if ~isempty(divTr)
            if tr < divTr && divTr-tr <= divSplitLimit(1) && divTr-tr > divSplitOmit(1)
                IDX_time_divSplit{k,1}(timestamps(tRange1):timestamps(tRange2)) = 1;
            elseif tr == divTr
                IDX_time_divSplit{k,1}(timestamps(tRange1):timestamps(tRange2)) = 2;
            elseif tr > divTr && tr-divTr <= divSplitLimit(2) && tr-divTr > divSplitOmit(2)
                IDX_time_divSplit{k,1}(timestamps(tRange1):timestamps(tRange2)) = 3;
            end
        end

        %IA-Policy Split Trials - Using CaitTransition
        IAFirst = find(IDX_set{2,5}{k,2}==1,1,'first');
        IALast = find(IDX_set{2,5}{k,2}==1,1,'last');
        if ~isempty(IAFirst) 
            if IDX_set{2,5}{k,2}(tr) == 1 && tr >= IAFirst && tr-IAFirst < IASplitLimit(1) && IALast-tr > IASplitMinDist %Must be IA-Policy; Must be on or after IAFirst, Must be within Limit(1) from IAFirst, must be further than MinDist from IALast
                IDX_time_IASplit{k,1}(timestamps(tRange1):timestamps(tRange2)) = 1;
            elseif IDX_set{2,5}{k,2}(tr) == 1 && tr <= IALast && IALast-tr < IASplitLimit(2) && tr-IAFirst > IASplitMinDist %Must be IA-Policy; Must be on or before IALast, Must be within Limit(1) from IALast, must be further than MinDist from IAFirst
                IDX_time_IASplit{k,1}(timestamps(tRange1):timestamps(tRange2)) = 2;
            elseif tr == divTr
                IDX_time_IASplit{k,1}(timestamps(tRange1):timestamps(tRange2)) = 3;
            end
        end

        

        %IDX_time_caitComb - Non-Phase = 0; IA-Policy = 1; Divergent(S+E) =2; Exploratory = 3; RS-Policy = 4;
        if IDX_set{2,16}{k,1}(tr)==1
            IDX_time_caitCombxPair{k,1}(timestamps(tRange1):timestamps(tRange2)) = IDX_set{2,4}{k,2}(tr);
        elseif IDX_set{2,16}{k,1}(tr)==2
            IDX_time_caitCombxPair{k,1}(timestamps(tRange1):timestamps(tRange2)) = IDX_set{2,4}{k,2}(tr)+4;
        end
            

    end

end

clear IDX_ALL_CAT IDX_ALL_INFO
IA=1; RS=2;

%Concatenate all matrix sets
IDX_ALL_INFO{1,1} = 'Title'; IDX_ALL_INFO{1,2} = 'IDX Names'; IDX_ALL_INFO{1,3} = 'IDX Values';
idx_cnt = 0;

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'Task Phase';  IDX_ALL_INFO{2,2}{idx_cnt} = {'Baseline','IA','RS'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1:3];
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_task{1:end});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'Trial Section'; IDX_ALL_INFO{2,2}{idx_cnt} = {'Pre-Dec','Post-Out','ITI'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1:3]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_phase{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'CE'; IDX_ALL_INFO{2,2}{idx_cnt} = {'Correct','Error'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1:2]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_CE{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'RS: Pair'; IDX_ALL_INFO{2,2}{idx_cnt} = {'Non-Conflict','Conflict'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1:2]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_pair{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'CaitClass'; IDX_ALL_INFO{2,2}{idx_cnt} = {'IA-Policy','Div(S+E)','Exploratory','RS-Policy'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1:4]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_caitComb{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'SubjectID'; IDX_ALL_INFO{2,2}{idx_cnt} = timepoints(:,1)'; IDX_ALL_INFO{2,3}{idx_cnt} = [1:size(timepoints,1)]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_subject{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'Manual Selection'; IDX_ALL_INFO{2,2}{idx_cnt} = tSelectionNames; IDX_ALL_INFO{2,3}{idx_cnt} = [1:size(tSelection,2)]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_selectTrials{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'CaitTransition'; IDX_ALL_INFO{2,2}{idx_cnt} = {'IA-Policy','IA-Div','Exploratory','RS-Acquire','RS-Policy'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1:5]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_caitTrans{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'DivergentSplit'; IDX_ALL_INFO{2,2}{idx_cnt} = {'Pre-Div','IA-Div','Post-Div'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1:3]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_divSplit{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'IAPolicySplit'; IDX_ALL_INFO{2,2}{idx_cnt} = {'FirstIA','FinalIA','IA-Div'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1:3]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_IASplit{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'All'; IDX_ALL_INFO{2,2}{idx_cnt} = {'All'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_All{1:end,1});

idx_cnt = idx_cnt+1;
IDX_ALL_INFO{2,1}{idx_cnt} = 'Trial'; IDX_ALL_INFO{2,2}{idx_cnt} = {'1...MaxTrials'}; IDX_ALL_INFO{2,3}{idx_cnt} = [1]; 
IDX_ALL_CAT(idx_cnt,:) = cat(2,IDX_time_trial{1:end,1});


clear IDX_TITLES
n=0;
for i=1:size(IDX_ALL_INFO{2,1},2)
    for j=1:size(IDX_ALL_INFO{2,2}{i},2)
        n=n+1;
        IDX_TITLES(n)=sprintf("%s || %s",IDX_ALL_INFO{2,1}{i},IDX_ALL_INFO{2,2}{i}{j});
    end
end

clear tRange1 tRange2 k n timestamps IA RS idx_cnt tr i j
clear IDX_time_task IDX_time_trial IDX_time_phase IDX_time_CE IDX_time_pCE IDX_time_pair IDX_time_ppair IDX_time_caitComb IDX_time_caitCombxPair
%%

%Step 5.1

clear stateSummary stateNames

for cstate=1:max(stateArray(1,:))
    n = 1;
    for IDX=1:size(IDX_ALL_CAT,1)
        clear temp
        for value=1:size(IDX_ALL_INFO{2,3}{IDX},2)
            temp(value) = round(sum(stateArray(1,:)==cstate & IDX_ALL_CAT(IDX,:)==IDX_ALL_INFO{2,3}{IDX}(value))/sum(IDX_ALL_CAT(IDX,:)==IDX_ALL_INFO{2,3}{IDX}(value))*100,0);
        end
        stateSummary(cstate,n:n+size(IDX_ALL_INFO{2,3}{IDX},2)-1) = round(temp/sum(temp)*100,0);
        n = n+size(IDX_ALL_INFO{2,3}{IDX},2);

    end

    stateNames(cstate)="";
    for i=1:size(CText,2)-2
        if ~isempty(CText{cstate,i})
            stateNames(cstate) = stateNames(cstate) + " | " + CText{cstate,i};
        end
    end
    if  stateNames(cstate)==""
         stateNames(cstate)="| No Synch";
    end

end

stateTable = array2table(stateSummary,"RowNames",stateNames,"VariableNames",IDX_TITLES);



filename = 'stateTable.xlsx';
writetable(stateTable,filename,'Sheet',1,'Range','D1')

clear value i j filename IDX temp n cstate

%% Step 5.2 - Build the synchTable
clear stateSummary synchSummary stateNames IDX_TITLES stateTable synchTable

averageSubjects = true;
useInclusiveSA = true; %Use the inclusive 'stateArrayFull' or the exclusive 'stateArray' to calculate

selection = 'divSplit'; %'default' or 'manual' or 'divSplit'or 'caitTrans' or 'IASplit' - manual utilizes the IDX_time_selectTrials instead of standard IDX combinations, divSplit is a unique case that examines pre and post divergent ONLY, caitTrans is strictly looking at classifications, not splitting further

%False + False - raw % of timepoints with state; True + False - %timepoints w/ state as a fraction of total IDX pieces; False + True - % timepoints w/ state as a fraction of all states
normalizeToSplit = false; %Calculates outputs as a fraction of a state's frequency / IDX within a split grouping
normalizeToIDX = false; %Calculates output as each state's fraction for each IDX combination

%Sets the IDX_SELECTION and IDX_RANGES based on the selection

%stateSet = [62 59 65 71];
stateSet = [41 31 38 30];

switch selection
    case 'default'
        IDX_SELECTION = [1 2 8 4 3]; %Version 6 is Subject(6) x Task Phase(1) x CaitTransition(8) x Pair(4) x CE(3)
        IDX_RANGES = {[3],[1:2],[1:5],[1:2],[1:2]};
    case 'manual'
        IDX_SELECTION = [1 2 7];
        IDX_RANGES = {[3],[1:2],IDX_ALL_INFO{2,3}{1,7}};
    case 'divSplit'
        IDX_SELECTION = [1 2 9 4];
        IDX_RANGES = {[3],[1],[3],[1,2]};
    case 'caitTrans'
        IDX_SELECTION = [1 2 8];
        IDX_RANGES = {[3],[1:2],[1:5]};
    case 'IASplit'
        IDX_SELECTION = [1 2 10 3];
        IDX_RANGES = {[3],[1:2],[1,2],[1:2]};
end


if averageSubjects == true
    numIter = max(IDX_ALL_INFO{2,3}{1,6});
    IDX_SELECTION = [6 IDX_SELECTION]; %Adds the Animal ID as the first SELECTION paramater
    IDX_RANGES = [{0} IDX_RANGES];
else
    numIter = 1;
end

for k=1:numIter
    fprintf('Performing aggregation %d/%d\n',k,numIter)

    if averageSubjects == true
        IDX_RANGES{1} = k;
    end

    OMIT_CE = [];
    %OMIT_CE = [2 1;8 1;8 5;9 1;9 2]; %Column 1: IDX_Selection, Column2: IDX_Range
    splitON = 2; %which index to split the table into groups

    T=combinations(IDX_RANGES{:});

    if strcmp(selection,'default') || strcmp(selection,'divSplit')
        for i=1:size(T,1)
            for j=1:size(OMIT_CE,1)
                idx = find(IDX_SELECTION==OMIT_CE(j,1));
                if T{i,idx}==OMIT_CE(j,2)
                    idx = find(IDX_SELECTION==3);
                    T{i,idx} = 99;
                end
            end
        end
    end
    T=unique(T); T=standardizeMissing(T,99); %Remove duplicate rows and then make 99s NaN - note: unique does not work with NaN values, thus this step is required
T3=T;
    for i=1:size(T,1)
        values = rmmissing(T{i,:});
        IDX_TITLES(i) = "";
        for j=1:size(values,2)
            if j==size(values,2)
                IDX_TITLES(i) = IDX_TITLES(i) + IDX_ALL_INFO{2,2}{IDX_SELECTION(j)}{values(j)};
            else
                IDX_TITLES(i) = IDX_TITLES(i) + IDX_ALL_INFO{2,2}{IDX_SELECTION(j)}{values(j)} + " | ";
            end
        end
    end

    for cstate=1:max(stateArray(1,:))
        m=1; n=0;
        clear temp temp2
        for i=1:size(T,1)
            n=n+1;
            values = rmmissing(T{i,:});
            combIDX = ones(1,size(IDX_ALL_CAT,2));
            for j=1:size(values,2)
                combIDX = combIDX & IDX_ALL_CAT(IDX_SELECTION(j),:)==IDX_ALL_INFO{2,3}{IDX_SELECTION(j)}(values(j));
            end

            if normalizeToIDX==true
                if useInclusiveSA==true
                    temp(n) = sum(stateArrayFull(cstate,:) & combIDX); %Total # of state occurances within IDX combination
                else
                    temp(n) = sum(stateArray(1,:)==cstate & combIDX); %Total # of state occurances within IDX combination
                end
            else
                if useInclusiveSA==true
                    temp(n) = sum(stateArrayFull(cstate,:) & combIDX)/sum(combIDX); %Percentage of state occurances within IDX combination
                else
                    temp(n) = sum(stateArray(1,:)==cstate & combIDX)/sum(combIDX); %Percentage of state occurances within IDX combination
                end

            end
            if useInclusiveSA==true
                temp2(n) = mean(stateArray(2,(stateArrayFull(cstate,:) & combIDX)));
            else
                temp2(n) = mean(stateArray(2,(stateArray(1,:)==cstate & combIDX)));
            end


            if i == size(T,1) || T{i,splitON}~=T{i+1,splitON}
                if normalizeToSplit==true
                    stateSummary(cstate,m:m+n-1) = round(temp/sum(temp,'omitnan'),4);
                else
                    stateSummary(cstate,m:m+n-1) = temp;
                end
                synchSummary(cstate,m:m+n-1) = temp2;
                m=m+n; n=0;
                clear temp
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Build the stateNames array%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        stateNames(cstate)="";
        for i=1:size(CText,2)-2
            if ~isempty(CText{cstate,i})
                stateNames(cstate) = stateNames(cstate) + " | " + CText{cstate,i};
            end
        end
        if  stateNames(cstate)==""
            stateNames(cstate)="| No Synch";
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    clear temp stateTableIDX
    if normalizeToIDX==true
        stateTableIDX = cell(2,size(stateSummary,2));
        for i=1:size(stateSummary,2)

            stateSummary(:,i) = stateSummary(:,i)/sum(stateSummary(:,i),1,"omitnan");

            [~,I] = sort(stateSummary(:,i),'descend');
            stateTableIDX{1,i} = IDX_TITLES(i);
            stateTableIDX{2,i} = array2table(stateSummary(I,i),"RowNames",stateNames(I),"VariableNames",IDX_TITLES(i));
        end
    end

    stateTable{k} = array2table(stateSummary,"RowNames",stateNames,"VariableNames",IDX_TITLES);
    synchTable{k} = array2table(synchSummary,"RowNames",stateNames,"VariableNames",IDX_TITLES);

    filename = 'stateTableDetailed.xlsx';
    writetable(stateTable{k},filename,'Sheet',k,'Range','D1')

    clear cstate i j m n o temp filename value1 value2 value3 value4 value5 T
end

clear stateTableMean stateTableSE
if averageSubjects == true
    temp = stateTable{1}{:,:};
    for k=2:9
        temp = cat(3,temp,stateTable{k}{:,:});
    end
    stateTableMean=array2table(mean(temp,3,'omitnan'),"RowNames",stateNames,"VariableNames",IDX_TITLES);
    stateTableSE=array2table(sem(temp,3),"RowNames",stateNames,"VariableNames",IDX_TITLES);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Step 6.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear grab
%1: Pre-Div NC, 2: Pre-Div C, 3: Post-Div NC, 4: Post-Div Conf
tblColumns = []; %Leave empty to grab all columns
%(AcrossHemisphere): 4=PVxPV; 16=+sHOOP; 20=+sHIP; 28=OPP; 42=xHOOP; 50=xHIP
%(WithinHemisphere): 4=PVxPV; 17=+sHOOP; 19=+shIP; 26=OPP; 43=xHOOP; 49=xHIP
n=0;
%stateSet = [62 59 65 71];
%stateSet = [99 80 89 94];
if isempty(tblColumns)
    tblColumns = [1:size(stateTable{1,1},2)];
end
for cstate = stateSet
    n = n+1;
    for section=tblColumns
        grab{section,1} = section;
        for k=1:size(stateTable,2)
            grab{section,2}(n,k) = stateTable{1,k}{cstate,section};
        end
    end
end

%grab = combCells(grab(:,2)','double','pre');



%% Polar plots based on state constraints
close all
plotThreshold = true;
plotPercent = true;
performWatson = true;

phaseBinIncrement = 15; %In degrees - How many degrees should each bin be incremented (Only for non-histogram plot as this allows for overlap)
phaseBinWidth = 60; %In degrees - how wide should the range of phases be for calculating inclusion (e.g. -45:45 is a length of 90)
histogramON = false; %Set to true for a traditional histogram, false for the polarplot version (allows for overlap between bins) 

oscChs = {'PV'}; %wH, xH, PV, or MD, wH-PD, xH-PD, wH-L, wH-R, xH-L, xH-R
constraints = {[3]}; %16, PVxPVExcl, 15 SHIP, 12 SHOOP, 42, xHIP, 35, xHOOP, 22 xHOPP
IDX_SELECTION = [1];
IDX_RANGES = {[3]};
OMIT_CE = []; %Column 1: IDX_Selection, Column2: IDX_Range



oscdata = cell.empty(0,5);
oscdata{1} = CText;
oscdata{2} = stateArrayFull;
if plotThreshold == true
    oscdata{3} = phasesRawAllThresh;
else
    oscdata{3} = phasesRawAll;
end
oscdata{4} = syncIn;
oscdata{5} = IDX_ALL_INFO;
oscdata{6} = IDX_ALL_CAT;


%output = plotoscphases(oscdata,oscChs,constraints,IDX_SELECTION,IDX_RANGES,OMIT_CE,plotPercent,plotThreshold,phaseBinIncrement,phaseBinWidth,histogramON,performWatson);
output2 = plotoscsynch(oscdata,oscChs,constraints,IDX_SELECTION,IDX_RANGES,OMIT_CE,plotThreshold,false,false,false);

% oscChSuperSet={{'wH-L','wH-R'},{'xH-L','xH-R'},{'MD'}};
% oscChSuperSetLabels = {'Ipsilateral PVIxPFC→MD Entrainment','Contralateral PVIxPFC→MD Entrainment','PFC→MD Interhemispheric Synchrony'};
% output3 = plotoscphaseheatmap(oscdata,oscChSuperSet,oscChSuperSetLabels,constraints,IDX_SELECTION,IDX_RANGES,OMIT_CE,plotPercent,plotThreshold,phaseBinIncrement,phaseBinWidth,histogramON,performWatson);

distFig


%% Loop to identify states that SIGNIFICANTLY enhance/suppress synchrony between two populations
close all
plotThreshold = false;
plotPercent = false;
performWatson = true;
averageSubjects = false;

phaseBinIncrement = 30; %In degrees - How many degrees should each bin be incremented (Only for non-histogram plot as this allows for overlap)
phaseBinWidth = 60; %In degrees - how wide should the range of phases be for calculating inclusion (e.g. -45:45 is a length of 90)
histogramON = true; %Set to true for a traditional histogram, false for the polatplot version (allows for overlap between bins)

oscdata = cell.empty(0,5);
oscdata{1} = CText;
oscdata{2} = stateArrayFull;
if plotThreshold == true
    oscdata{3} = phasesRawAllThresh;
else
    oscdata{3} = phasesRawAll;
end
oscdata{4} = syncIn; %Sets the IN-PHASE synchrony as the synchrony data
oscdata{5} = IDX_ALL_INFO;
oscdata{6} = IDX_ALL_CAT;

OMIT_CE = []; %Column 1: IDX_Selection, Column2: IDX_Range
oscChs = {'MD'}; %wH, xH, PV, or MD, wH-PD, xH-PD, wH-L, wH-R, xH-L, xH-R

CTRLSTATE = 4; CTRLCLMN = 1;  %Row (CTRLSTATE) and Column (CTRLCLMN) in CText that identifies the control state to use in the loop - note: should only select states with a single active column

IDX_SELECTION = [1 9];
IDX_RANGES = {[3],[3]};

if averageSubjects == true
    numIter = max(IDX_ALL_INFO{2,3}{1,6});
    IDX_SELECTION = [6 IDX_SELECTION]; %Adds the Animal ID as the first SELECTION paramater
    IDX_RANGES = [{0} IDX_RANGES];
else
    numIter = 1;
end

sigStateOUT = cell.empty(0,2);
sigStateOUT{1} = string.empty(0,length(find(strcmp(CText(:,1),CText{CTRLSTATE,CTRLCLMN}))));
sigStateOUT{2} = double.empty(0,length(find(strcmp(CText(:,1),CText{CTRLSTATE,CTRLCLMN})))*3);

for k=1:numIter
    if averageSubjects == true
        IDX_RANGES{1} = k;
    end

    o=0; n=0; m=0; l=0; clear sigStates sigStateIncrease sigStateDecrease sigStateNoChange

    for cstate=1:size(CText,1)
        if strcmp(CText{cstate,CTRLCLMN},CText{CTRLSTATE,CTRLCLMN}) && cstate ~= CTRLSTATE

            constraints = {[CTRLSTATE],[cstate]}; %2 10 - RH, 2 12 -LH
            output2 = plotoscsynch(oscdata,oscChs,constraints,IDX_SELECTION,IDX_RANGES,OMIT_CE,plotThreshold,false,plotPercent,true);
            if o == 0 %First state should be the control state
                sigStateOUT{1}(1) = sprintf('State %d',CTRLSTATE);
                sigStateOUT{2}(k,1) = mean(output2{1,2});
                sigStateOUT{2}(k,2) = std(output2{1,2})/sqrt(length(output2{1,2}));
                sigStateOUT{2}(k,3) = length(output2{1,2});
                o = o+1;
            end
            sigStateOUT{1}(o+1) = sprintf('State %d',cstate);
            sigStateOUT{2}(k,o*3 + 1) = mean(output2{2,2});
            sigStateOUT{2}(k,o*3 + 2) = std(output2{2,2})/sqrt(length(output2{2,2}));
            sigStateOUT{2}(k,o*3 + 3) = length(output2{2,2});
            o = o+1;

            [h,p] = ttest2(output2{1,2},output2{2,2});
            if p<=0.05 && mean(output2{1,2}) > mean(output2{2,2})
                n=n+1;
                sigStateDecrease(n,1) = cstate;
                sigStateDecrease(n,2) = o;
                sigStateDecrease(n,3) = p;
                sigStateDecrease(n,4) = mean(output2{2,2})/mean(output2{1,2})*100;
            elseif p<=0.05 && mean(output2{1,2}) < mean(output2{2,2})
                m=m+1;
                sigStateIncrease(m,1) = cstate;
                sigStateIncrease(m,2) = o;
                sigStateIncrease(m,3) = p;
                sigStateIncrease(m,4) = mean(output2{2,2})/mean(output2{1,2})*100;
            else
                l=l+1;
                sigStateNoChange(l,1) = cstate;
                sigStateNoChange(l,2) = o;
                sigStateNoChange(l,3) = p;
                sigStateNoChange(l,4) = mean(output2{2,2})/mean(output2{1,2})*100;
            end
            close all
        end
    end
end
