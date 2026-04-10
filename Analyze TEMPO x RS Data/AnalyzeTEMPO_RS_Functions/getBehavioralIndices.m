function [IDX_set,metadata] = getBehavioralIndices(timepoints)

timepoints = [timepoints(:,1) timepoints(:,1) timepoints(:,2:end)];
AnalyzeRSStrategies
timepoints = [timepoints(:,1) timepoints(:,3:end)];
clear NRB_NC NRB_OC Sidx

%Based on over 25 mice, less than 8% mice have 2 errors in the last 6 trials, however this jumps
%up to 20% at last 7 trials. Thus the last 5 or 6 trials should be
%considered the 'final' streak. Conversely, 32% mice have 0 errors in the
%last 6 trials, this drops to 16% in the last 7 trials, 12% in last 8 trials and 0% in last 9 trials.
%Therefore, the 'peri-streak' seems to end between the last 9 to 7 trials.
%Since this range can encompass critically different activity than
%post-strk, peri-streak should be determined flexibly: the 2 trials before
%and the 2 trials after the 2nd to last error. In the case of 1 error, the
%first 3 trials will be used.

%Behavioral trial ranges of interest, trials from end of behavior
rng_prestrk = [1 3]; %From Trial Start 1 = 1st trial
rng_poststrk = [7 0];%Distance from trial end 0 = last trial %[4 0] 7/3/23, temp change to old def of peristreak
rng_strkstart = 10;
rng_peristrk = [11 7]; %Distance from last error trial, [left right] shift; original peristreak from final trial: [11 7]
rng_peristrat = [3 0]; %Distance from last trial in a strategy block (either IDX_4BlockClassification or IDX_4BlockStages)
rng_IAEnd = [4 0];
rng_IAStart = [1 3];

%Limit to subset of correct or error trials
MF_CE_SubsetON = false;
MF_CE_range = 1; %# of errors/correct trials
MF_CE_Order = 'first'; %Either 'first' or 'last'

%Caitriona Classification variables
MF_CaitClass_CmbnDiv = true; %Set to true if you want to combine divergent trial types S->S and S->E
MF_CaitClass_DivExp2Stble = true; %If true, the trial preceding stable policies (except the first trial of IA/RS) will be set to divergent. - A divergence from 'exploratory' behavior

CaitClassSUBSET = 'During'; %'Before', 'During', 'After'; Sets aCaitclass to a subset in relation to the Behavioral Transition period

MF_CaitClass_subsetON = false; %Whether to limit the classifications to a subset of trials
MF_CaitClass_Range = 1; %Range of trials to include in subset
MF_CaitClass_Order = 'first'; %Either 'first' # of trials, or 'last' # of trials

MF_CaitStrat_LimitBlock = false; %Includes only the 'First' or 'Final' policy blocks.
MF_CaitStrat_LimitBlock_Order = 'last'; %first or last

MF_CaitStrat_subsetON = false; %Same as above, but for the stable policies
MF_CaitStrat_Range = 3;
MF_CaitStrat_Order = 'last'; %first, last, before, or after

%Metadata Structure
metadata = struct('rng_prestrk',rng_prestrk,'rng_poststrk',rng_poststrk,'rng_strkstart',rng_strkstart,'rng_peristrk',rng_peristrk,'rng_peristrat',rng_peristrat,'rng_IAEnd',rng_IAEnd,'rng_IAStart',rng_IAStart,...
    'MF_CE_SubsetON',MF_CE_SubsetON,'MF_CE_range',MF_CE_range,'MF_CE_Order',MF_CE_Order,...
    'MF_CaitClass_CmbnDiv',MF_CaitClass_CmbnDiv,'MF_CaitClass_DivExp2Stble',MF_CaitClass_DivExp2Stble,'CaitClassSUBSET',CaitClassSUBSET,'MF_CaitClass_subsetON',MF_CaitClass_subsetON,'MF_CaitClass_Range',MF_CaitClass_Range,'MF_CaitClass_Order',MF_CaitClass_Order,'MF_CaitStrat_LimitBlock',MF_CaitStrat_LimitBlock,'MF_CaitStrat_LimitBlock_Order',MF_CaitStrat_LimitBlock_Order,'MF_CaitStrat_subsetON',MF_CaitStrat_subsetON,'MF_CaitStrat_Range',MF_CaitStrat_Range,'MF_CaitStrat_Order',MF_CaitStrat_Order);


%% Generate Dig Delay data
clear digdelays digdelays_zscore
for k=1:size(timepoints,1)
    %IA Dig Delays
    for tr=1:size(timepoints{k,2},1)
        digdelays{k,1}(tr,1) = timepoints{k,2}{tr,2}-timepoints{k,2}{tr,1};
    end
    %RS Dig Delays
    for tr=1:size(timepoints{k,3},1)
        digdelays{k,2}(tr,1) = timepoints{k,3}{tr,2}-timepoints{k,3}{tr,1};
    end

    digdelays{k,1}(:,2) = [digdelays{k,1}(2:end);nan]; %2nd Column is shifted for next trial's dig delay
    digdelays{k,2}(:,2) = [digdelays{k,2}(2:end);nan];

end



%Then get a z-scored (across all an animal's trials) version
clear temp
for k=1:size(digdelays,1)
    temp = [digdelays{k,1}(:,1);digdelays{k,2}(:,1)]; %Wanting to z-score across all trials, not differentially bw IA and RS
    temp = zscore(temp);
    digdelays_zscore{k,1} = temp(1:size(digdelays{k,1},1));
    digdelays_zscore{k,2} = temp(size(digdelays{k,1},1)+1:end);

    digdelays_zscore{k,1}(:,2) = [digdelays_zscore{k,1}(2:end);nan]; %2nd Column is shifted for next trial's dig delay
    digdelays_zscore{k,2}(:,2) = [digdelays_zscore{k,2}(2:end);nan];

end

%% Generate the IDX_set - Exhaustive List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear IDX_set
IDX_counter=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_AllTrials%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear IDX_AllTrials
for k=1:size(timepoints,1)

    IDX_AllTrials{k,1} = ones(size(timepoints{k,2},1) ,1);
    IDX_AllTrials{k,2} = ones(size(timepoints{k,3},1) ,1);

end

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'All Trials';
IDX_set{2,IDX_counter} = IDX_AllTrials;
IDX_set{3,IDX_counter} = {'All Trials'};
IDX_set{4,IDX_counter} = [1];
IDX_set{5,IDX_counter} = 'No grouping of trials - Plot everything';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_animal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IDX_animal lists each animal's ID # for all of it's trials in IA and RS
clear IDX_animal
for k=1:size(timepoints,1)
    for tr=1:size(timepoints{k,2},1)
        IDX_animal{k,1}(tr,1) = k;
    end
    for tr=1:size(timepoints{k,3},1)
        IDX_animal{k,2}(tr,1) = k;
    end
end

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'Animal';
IDX_set{2,IDX_counter} = IDX_animal;
IDX_set{3,IDX_counter} = {timepoints(:,1)'};
IDX_set{4,IDX_counter} = [1:size(timepoints,1)];
IDX_set{5,IDX_counter} = 'Examine each animal individually';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_CE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear IDX_CE
for k=1:size(timepoints,1)

    %Initial-Association
    for tr=1:size(timepoints{k,2},1)
        IDX_CE{k,1}(tr,1) = timepoints{k,2}{tr,6};
    end

    %Rule-Shift
    for tr=1:size(timepoints{k,3},1)
        IDX_CE{k,2}(tr,1) = timepoints{k,3}{tr,6};
    end


    if MF_CE_SubsetON==true
        idx1 = find(IDX_CE{k,1}==0,MF_CE_range,MF_CE_Order);
        idx2 = find(IDX_CE{k,1}==1,MF_CE_range,MF_CE_Order);
        IDX_CE{k,1}(:)=NaN;
        IDX_CE{k,1}(idx1) = 0;
        IDX_CE{k,1}(idx2) = 1;

        idx1 = find(IDX_CE{k,2}==0,MF_CE_range,MF_CE_Order);
        idx2 = find(IDX_CE{k,2}==1,MF_CE_range,MF_CE_Order);
        IDX_CE{k,2}(:)=NaN;
        IDX_CE{k,2}(idx1) = 0;
        IDX_CE{k,2}(idx2) = 1;
    end


end

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'Outcome';
IDX_set{2,IDX_counter} = IDX_CE;
IDX_set{3,IDX_counter} = {'Correct','Error'};
IDX_set{4,IDX_counter} = [0:1];
IDX_set{5,IDX_counter} = 'Split between outcome type: Correct vs Error';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_pair%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear IDX_pair %Reinforced Pair = 1; New Pair(i.e. requires new strat) = 2

for k=1:size(timepoints,1)

    for tr=1:size(timepoints{k,2},1) %No IDX in IA - all NaN - SHOULD back-analyze and set the conflict/nonconflict pairs for IA as well
        IDX_pair{k,1}(tr,1) = NaN;
    end

    for tr=1:size(timepoints{k,3},1)
        if  timepoints{k,3}{tr,6}==0 && timepoints{k,3}{tr,7}==1 %Correct and OS Required (ReinfPair)
            IDX_pair{k,2}(tr,1) = 1;
        elseif timepoints{k,3}{tr,6}==1 && timepoints{k,3}{tr,7}==1 %Error and OS Required (NewPair) - Random Error
            IDX_pair{k,2}(tr,1) = 2;
        elseif timepoints{k,3}{tr,6}==0 && timepoints{k,3}{tr,7}==0 %Correct and NS Required (NewPair)
            IDX_pair{k,2}(tr,1) = 2;
        elseif timepoints{k,3}{tr,6}==1 && timepoints{k,3}{tr,7}==0 %Error and NS Required (ReinfPair) - Perseverative Error
            IDX_pair{k,2}(tr,1) = 1;
        end
    end
end


IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'Pair Presentation';
IDX_set{2,IDX_counter} = IDX_pair;
IDX_set{3,IDX_counter} = {'Non-Conflict','Conflict'};
IDX_set{4,IDX_counter} = [1:2];
IDX_set{5,IDX_counter} = 'Split between pair presentations: Non-Conflict and Conflict';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_digdelay%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Uses the calculated dig delay (e.g. time from trial start to dig start) to
%classify trials - based on z-scores within animal
clear IDX_digdelay
for k=1:size(timepoints,1)

    IDX_digdelay{k,1} = zeros(size(timepoints{k,2},1) ,1);
    IDX_digdelay{k,2} = zeros(size(timepoints{k,3},1) ,1);

    %First IA
    for tr=1:size(timepoints{k,2},1)
        if digdelays_zscore{k,1}(tr,1) > 0.5 %Slow Decision
            IDX_digdelay{k,1}(tr,1) = 2;
        elseif digdelays_zscore{k,1}(tr,1) < -0.5 %Fast Decision
            IDX_digdelay{k,1}(tr,1) = 0;
        else %Normal Decision
            IDX_digdelay{k,1}(tr,1) = 1;
        end
    end
    %Next RS
    for tr=1:size(timepoints{k,3},1)
        if digdelays_zscore{k,2}(tr,1) > 0.5 %Slow Decision
            IDX_digdelay{k,2}(tr,1) = 2;
        elseif digdelays_zscore{k,2}(tr,1) < -0.5 %Fast Decision
            IDX_digdelay{k,2}(tr,1) = 0;
        else %Normal Decision
            IDX_digdelay{k,2}(tr,1) = 1;
        end
    end
end


IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'Dig Delay';
IDX_set{2,IDX_counter} = IDX_digdelay;
IDX_set{3,IDX_counter} = {'Fast','Normal','Slow'};
IDX_set{4,IDX_counter} = [0:2];
IDX_set{5,IDX_counter} = 'Separate fast, normal, and slow delays from trial start to dig';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_chronological%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IDX_chronological are chronological parcelation of trials based on criterion
%defined above
clear IDX_chronological %1 = Early Behavior (Within rng_prestrk) 2 = Late Behavior (within rng_strk)
for k=1:size(timepoints,1)
    %IA
    IDX_chronological{k,1} = zeros(size(timepoints{k,2},1) ,1);

    if size(timepoints{k,2},1) > rng_IAStart(2) + rng_strkstart
        IDX_chronological{k,1}(rng_IAStart(1):rng_IAStart(2)) = 1;
        IDX_chronological{k,1}(end-rng_IAEnd(1)+1:end-rng_IAEnd(2)) = 3;
    elseif size(timepoints{k,2},1) >= rng_IAStart(1) + rng_strkstart + 2
        IDX_chronological{k,1}(rng_IAStart(1):end-rng_strkstart) = 1;
        IDX_chronological{k,1}(end-rng_IAEnd(1)+1:end-rng_IAEnd(2)) = 3;
    elseif size(timepoints{k,2},1) > rng_strkstart
        IDX_chronological{k,1}(1:3) = 1;
        IDX_chronological{k,1}(end-rng_IAEnd(1)+1:end-rng_IAEnd(2)) = 3;
    elseif size(timepoints{k,2},1) > rng_IAEnd(1)
        IDX_chronological{k,1}(1:3) = 1;
        IDX_chronological{k,1}(end-rng_IAEnd(1)+1:end-rng_IAEnd(2)) = 3;
    elseif size(timepoints{k,2},1) >= rng_IAEnd(2)
        IDX_chronological{k,1}(1:3) = 1;
        IDX_chronological{k,1}(1:end-rng_IAEnd(2)) = 3;
    elseif size(timepoints{k,2},1) < rng_IAEnd(2)
        %No trials in range, leave all as 0
    end

    IDX_chronological{k,1}(IDX_chronological{k,1}==0) = 2;

    %RS
    IDX_chronological{k,2} = zeros(size(timepoints{k,3},1) ,1);

    if size(timepoints{k,3},1) > rng_prestrk(2) + rng_strkstart
        IDX_chronological{k,2}(rng_prestrk(1):rng_prestrk(2)) = 1;
        IDX_chronological{k,2}(end-rng_poststrk(1)+1:end-rng_poststrk(2)) = 3;
    elseif size(timepoints{k,3},1) >= rng_prestrk(1) + rng_strkstart + 2
        IDX_chronological{k,2}(rng_prestrk(1):end-rng_strkstart) = 1;
        IDX_chronological{k,2}(end-rng_poststrk(1)+1:end-rng_poststrk(2)) = 3;
    elseif size(timepoints{k,3},1) > rng_strkstart
        IDX_chronological{k,2}(1:3) = 1;
        IDX_chronological{k,2}(end-rng_poststrk(1)+1:end-rng_poststrk(2)) = 3;
    elseif size(timepoints{k,3},1) > rng_poststrk(1)
        IDX_chronological{k,2}(1:3) = 1;
        IDX_chronological{k,2}(end-rng_poststrk(1)+1:end-rng_poststrk(2)) = 3;
    elseif size(timepoints{k,3},1) >= rng_poststrk(2)
        IDX_chronological{k,2}(1:3) = 1;
        IDX_chronological{k,2}(1:end-rng_poststrk(2)) = 3;
    elseif size(timepoints{k,3},1) < rng_poststrk(2)
        %No trials in range, leave all as 0
    end

    IDX_chronological{k,2}(IDX_chronological{k,2}==0) = 2;

end

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'Chronological Ranges';
IDX_set{2,IDX_counter} = IDX_chronological;
IDX_set{3,IDX_counter} = {'Early','Mid','Late'};
IDX_set{4,IDX_counter} = [1:3];
IDX_set{5,IDX_counter} = 'Early, mid, versus late behavior';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_peristreak%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Peri-streak based on end of behavior, not final error
clear IDX_peristreak %1 = Early Behavior (Within rng_prestrk) 2 = Late Behavior (within rng_strk)
for k=1:size(timepoints,1)
    %IA
    IDX_peristreak{k,1} = zeros(size(timepoints{k,2},1) ,1);
    IDX_peristreak{k,2} = zeros(size(timepoints{k,3},1) ,1);

    %Peri-streak period
    if size(timepoints{k,2},1) > rng_prestrk(2) + rng_strkstart
        IDX_peristreak{k,1}(end-rng_peristrk(1)+1:end-rng_peristrk(2)) = 1;
    elseif size(timepoints{k,2},1) >= rng_prestrk(1) + rng_strkstart + 2
        IDX_peristreak{k,1}(end-rng_peristrk(1)+1:end-rng_peristrk(2)) = 1;
    elseif size(timepoints{k,2},1) > rng_peristrk(1)
        IDX_peristreak{k,1}(end-rng_peristrk(1)+1:end-rng_peristrk(2)) = 1;
    elseif size(timepoints{k,2},1) > rng_strkstart
        IDX_peristreak{k,1}(end-rng_strkstart+1:end-rng_peristrk(2)) = 1;
    elseif size(timepoints{k,2},1) >= rng_peristrk(2)
        IDX_peristreak{k,1}(1:end-rng_peristrk(2)) = 1;
    elseif size(timepoints{k,2},1) < rng_peristrk(2)
        %No trials in range, leave all as 0
    end

    %RS

    IDX_peristreak{k,2} = zeros(size(timepoints{k,3},1) ,1);

    %Peri-streak period
    if size(timepoints{k,3},1) > rng_prestrk(2) + rng_strkstart
        IDX_peristreak{k,2}(end-rng_peristrk(1)+1:end-rng_peristrk(2)) = 1;
    elseif size(timepoints{k,3},1) >= rng_prestrk(1) + rng_strkstart + 2
        IDX_peristreak{k,2}(end-rng_peristrk(1)+1:end-rng_peristrk(2)) = 1;
    elseif size(timepoints{k,3},1) > rng_peristrk(1)
        IDX_peristreak{k,2}(end-rng_peristrk(1)+1:end-rng_peristrk(2)) = 1;

    elseif size(timepoints{k,3},1) > rng_strkstart
        IDX_peristreak{k,2}(end-rng_strkstart+1:end-rng_peristrk(2)) = 1;
    elseif size(timepoints{k,3},1) >= rng_peristrk(2)
        IDX_peristreak{k,2}(1:end-rng_peristrk(2)) = 1;
    elseif size(timepoints{k,3},1) < rng_peristrk(2)
        %No trials in range, leave all as 0
    end


    %idx = find(cell2mat(timepoints{k,3}(:,6))==1,3,"last");
    %IDX_peristreak{k,2}(idx(1):idx(end)) = 1;

end

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'Peri-Streak';
IDX_set{2,IDX_counter} = IDX_peristreak;
IDX_set{3,IDX_counter} = {'PeriStrk'};
IDX_set{4,IDX_counter} = [1];
IDX_set{5,IDX_counter} = 'Trials surrounding the start of the streak';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_4BlockClassification%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IDX_4BlockClassification uses AnalyzeRSStrategies to identify trials based on the
%supervised rolling 4-block trial classification method
clear IDX_4BlockClassification
for k=1:size(timepoints,1)
    for tr=1:size(timepoints{k,2},1)
        IDX_4BlockClassification{k,1}(tr,1) = NaN;
    end
    IDX_4BlockClassification{k,2} = combStrat{1,k};
end

%%%%%%%%%%%%%%%%%%IDX_4BlockStages%%%%%%%%%%%%%%%%%%
%IDX_4BlockStages is an expansion of IDX_4BlockClassification but generalizes strategies into
%three categories
clear IDX_4BlockStages %1 = Perseveration 2 = Flexible Learning 3 = New Strategy Implemented (defined as trials after final error)
for k=1:size(timepoints,1)

    for tr=1:size(timepoints{k,2},1) %No IDX in IA - all NaN
        IDX_4BlockStages{k,1}(tr,1) = NaN;
    end

    IDX_4BlockStages{k,2} = zeros(size(timepoints{k,3},1) ,1);

    IDX_4BlockStages{k,2}(IDX_4BlockClassification{k,2}==0) = 1; %Perseveration as defined by strat for now
    IDX_4BlockStages{k,2}(IDX_4BlockClassification{k,2}==1) = 99; %Regressive behavior excluded
    IDX_4BlockStages{k,2}(IDX_4BlockClassification{k,2}==3) = 99; %Reversal behavior excluded

    %     if sum(IDX_4BlockClassification{k,2}==1,1) > 0 %Animal fell into regressive behavior
    %             IDX_4BlockStages{k,2}(find(IDX_4BlockClassification{k,2}==1,2,'last'))= 2; % Last 2 perseveration included in flexible learning
    %     else
    %         if sum(IDX_4BlockClassification{k,2}==0,1) > 2
    %             IDX_4BlockStages{k,2}(find(IDX_4BlockClassification{k,2}==0,2,'last'))= 2; % Last 2 perseveration included in flexible learning
    %         elseif sum(IDX_4BlockClassification{k,2}==0,1) > 1
    %             IDX_4BlockStages{k,2}(find(IDX_4BlockClassification{k,2}==0,1,'last'))= 2; % Last 1 perseveration included in flexible learning
    %         end
    %     end

    if size(IDX_4BlockClassification{k,2},1) > 10
        if sum(IDX_4BlockClassification{k,2}==0,1) > 2
            IDX_4BlockStages{k,2}(find(IDX_4BlockClassification{k,2}==0,2,'last'))= 2; % Last 2 perseveration included in flexible learning
        elseif sum(IDX_4BlockClassification{k,2}==0,1) > 1
            IDX_4BlockStages{k,2}(find(IDX_4BlockClassification{k,2}==0,1,'last'))= 2; % Last 1 perseveration included in flexible learning
        end
        %IDX_4BlockStages{k,2}(find(IDX_4BlockClassification{k,2}==1,2,'last'))= 2; % Last 2 regressive included in flexible learning
    end

    %     if sum(IDX_4BlockClassification{k,2}==1,1) > 1
    %         IDX_4BlockStages{k,2}(find(IDX_4BlockClassification{k,2}==1,2,'last'))= 2; % Last 2 perseveration included in flexible learning
    %     end


    %%%%%%Applied Strategy%%%%%%
    %     if size(IDX_4BlockClassification{k,2},1) >= 10 %Last 10 trials is 'Applied Strategy'
    %         IDX_4BlockStages{k,2}(10:end) = 3;
    %     else
    %         IDX_4BlockStages{k,2}(1:end) = 3;
    %     end

    % IDX_4BlockStages{k,2}(find(IDX_4BlockClassification{k,2}==4,1,'first'):end) = 3;

    idx=find(cell2mat(timepoints{k,3}(:,6))==1,2,"last");
    IDX_4BlockStages{k,2}(idx(1)+1:end) = 3; %Define New Strategy Implemented as all trials after

    %      idx = find(cell2mat(timepoints{k,3}(:,6))==1 & ~IDX_4BlockClassification{k,2}==0,3,"last");%find(cell2mat(timepoints{k,3}(:,6))==1 & cell2mat(timepoints{k,3}(:,7))==1,1,"last"); %Find the last error trial
    %      if isempty(idx)
    %         IDX_4BlockStages{k,2}(find(IDX_4BlockClassification{k,2}==4,1,'first'):end) = 3;
    %      else
    %         IDX_4BlockStages{k,2}(idx(1):end) = 3; %Define New Strategy Implemented as all trials after
    %      end

    IDX_4BlockStages{k,2}(IDX_4BlockStages{k,2}==0)=2; %All remaining trials as 'Flexible Learning' period
end

clear IDX_4BlockBehTrans
for k=1:size(timepoints,1)

    for tr=1:size(timepoints{k,2},1) %No IDX in IA - all NaN
        IDX_4BlockBehTrans{k,1}(tr,1) = NaN;
    end

    IDX_4BlockBehTrans{k,2} = zeros(size(timepoints{k,3},1) ,1);

    clear idx
    idx = find(cell2mat(timepoints{k,3}(:,6))==1 & cell2mat(timepoints{k,3}(:,7))==1 & IDX_4BlockClassification{k,2}~=0,2,'last'); %Find 2nd to last conflict error not in perseverative blocks
    idx = sort([idx;find(cell2mat(timepoints{k,3}(:,6))==0 & cell2mat(timepoints{k,3}(:,7))==0,2,'first')]); %Find first 2 conflict correct trials

    % clear idx
    % idx(1) = find(IDX_Event{k,2}==5,1,'first');
    % idx(2) = find(IDX_Event{k,2}==5,1,'last');

    IDX_4BlockBehTrans{k,2}(idx(1):idx(end)) = 1;
    IDX_4BlockBehTrans{k,2}(idx(end)+1:end) = 2;

end

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = '4-Block Strategy Classification';
IDX_set{2,IDX_counter} = IDX_4BlockClassification;
IDX_set{3,IDX_counter} = {'Pers','Reg','NRB','Rev','Strk'};
IDX_set{4,IDX_counter} = [0:4];
IDX_set{5,IDX_counter} = 'Uses a rolling 4-block method to classify behavior into: Perseverative, Regressive, New Random Behavior, Reversal, and Final Streak';

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = '4-Block Learning Stages';
IDX_set{2,IDX_counter} = IDX_4BlockStages;
IDX_set{3,IDX_counter} = {'Perseveration','Flexible','Adapted Strat'};
IDX_set{4,IDX_counter} = [1:3];
IDX_set{5,IDX_counter} = 'Groups 4-Block classifications into three learning stages: Perseveration, Flexible Learning, and Adapted Strategy';

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = '4-Block Behavioral Transition';
IDX_set{2,IDX_counter} = IDX_4BlockBehTrans;
IDX_set{3,IDX_counter} = {'Pre-Transition', 'Behavioral Transition', 'Post-Transition'};
IDX_set{4,IDX_counter} = [0:2];
IDX_set{5,IDX_counter} = 'Groups 4-Block classifications into three learning phases: Behavioral Transition, Pre-Transition, and Post-Transition';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_CaitClass-IDX_CaitCombined2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The IDX_Cait set uses Caitriona's ML classification system to parcel
%trials into different classified categories - highly modifiable at the
%moment
clear IDX_CaitClass IDX_CaitPolicies IDX_CaitCombined1 IDX_CaitCombined2 
for k=1:size(timepoints,1)


    if size(timepoints,2) > 3
        tempIDX = timepoints{1,4}{1};
        IDX_CaitClass{k,1} = tempIDX{k,1}(1:size(timepoints{k,2},1));
        IDX_CaitClass{k,2} = tempIDX{k,2}(1:size(timepoints{k,3},1));
    else
        IDX_CaitClass{k,1} = ones(size(timepoints{k,2},1),1); %%IF THERE IS NOT A CAITCLASS - JUST MAKE ALL 1s
        IDX_CaitClass{k,2} = ones(size(timepoints{k,3},1),1); %%IF THERE IS NOT A CAITCLASS - JUST MAKE ALL 1s
    end

    if size(timepoints,2) > 3
        tempIDX = timepoints{1,4}{2};
        IDX_CaitPolicies{k,1} = tempIDX{k,1}(1:size(timepoints{k,2},1));
        IDX_CaitPolicies{k,2} = tempIDX{k,2}(1:size(timepoints{k,3},1));
    else
        IDX_CaitPolicies{k,1} = ones(size(timepoints{k,2},1),1); %%IF THERE IS NOT A CAITCLASS - JUST MAKE ALL 1s
        IDX_CaitPolicies{k,2} = ones(size(timepoints{k,3},1),1); %%IF THERE IS NOT A CAITCLASS - JUST MAKE ALL 1s
    end


    if MF_CaitClass_CmbnDiv == true
        IDX_CaitClass{k,2}(IDX_CaitClass{k,2}==3)=2; %Combine S->S and S->E Divergent but not Lapse
    end


    if MF_CaitClass_DivExp2Stble == true %Sets exploratory trials that are tr-1 of stable policies as divegent, i.e. they are diverging from exploratory to stable

        idx = find(IDX_CaitClass{k,2}==1); %All Stable policies
        idx = idx([1;flip(diff(flip(idx))~=-1)]==true); %Remove consecutive trials within a stable policy block (keep only the 1st trial)

        for tr=idx'-1 %For all preceding trials
            if tr > 0 && IDX_CaitClass{k,2}(tr)==0 %If trial exists (>0) and is an exploratory trial, mark as divergent
                if MF_CaitClass_CmbnDiv == true
                    IDX_CaitClass{k,2}(tr)=2;
                else
                    IDX_CaitClass{k,2}(tr)=5;
                end
            end
        end
    end

    %%%%%Only examine a subset of classifications (e.g. only 1st 3).
    if MF_CaitClass_subsetON == true
        subrange = MF_CaitClass_Range;
        for i=0:4
            idx = find(IDX_CaitClass{k,1}==i);
            if length(idx) > subrange
                if strcmp(MF_CaitClass_Order,'first')
                    IDX_CaitClass{k,1}(idx(subrange+1:end))=NaN;
                else
                    IDX_CaitClass{k,1}(idx(1:end-subrange))=NaN;
                end
            end

            idx = find(IDX_CaitClass{k,2}==i);
            if length(idx) > subrange
                if strcmp(MF_CaitClass_Order,'first')
                    IDX_CaitClass{k,2}(idx(subrange+1:end))=NaN;
                else
                    IDX_CaitClass{k,2}(idx(1:end-subrange))=NaN;
                end
            end
        end
    end

    %%%%Since CaitStrats occur in blocks, limit to the first or last block only
    if MF_CaitStrat_LimitBlock == true
        for i=0:6
            if strcmp(MF_CaitStrat_LimitBlock_Order,'first')

                if sum(IDX_CaitPolicies{k,1}==i)>0
                    mask = IDX_CaitPolicies{k,1}==i;
                    mask2 = IDX_CaitPolicies{k,1}==i;
                    idx1 = find(mask2==1);
                    idx2 = find(mask2==0);
                    mask2(idx2(find(idx2>idx1(1),1)):end) = 0;
                    IDX_CaitPolicies{k,1}(mask & ~mask2)=99;
                end

                if sum(IDX_CaitPolicies{k,2}==i)>0
                    mask = IDX_CaitPolicies{k,2}==i;
                    mask2 = IDX_CaitPolicies{k,2}==i;
                    idx1 = find(mask2==1);
                    idx2 = find(mask2==0);
                    mask2(idx2(find(idx2>idx1(1),1)):end) = 0;
                    IDX_CaitPolicies{k,2}(mask & ~mask2)=99;
                end
            end

            if strcmp(MF_CaitStrat_LimitBlock_Order,'last')

                if sum(IDX_CaitPolicies{k,1}==i)>0
                    mask = IDX_CaitPolicies{k,1}==i;
                    mask2 = IDX_CaitPolicies{k,1}==i;
                    idx1 = find(mask2==1);
                    idx2 = find(mask2==0);
                    mask2(1:idx2(find(idx2<idx1(end),1,'last'))) = 0;
                    IDX_CaitPolicies{k,1}(mask & ~mask2)=99;
                end

                if sum(IDX_CaitPolicies{k,2}==i)>0
                    mask = IDX_CaitPolicies{k,2}==i;
                    mask2 = IDX_CaitPolicies{k,2}==i;
                    idx1 = find(mask2==1);
                    idx2 = find(mask2==0);
                    mask2(1:idx2(find(idx2<idx1(end),1,'last'))) = 0;
                    IDX_CaitPolicies{k,2}(mask & ~mask2)=99;
                end
            end

        end
    end

    %%%%%%Only examine a subset of stable policies (e.g. only 1st 3)
    if MF_CaitStrat_subsetON == true
        subrange = MF_CaitStrat_Range;
        for i=0:6
            idx = find(IDX_CaitPolicies{k,1}==i);
            if length(idx) > subrange
                if strcmp(MF_CaitStrat_Order,'first')
                    IDX_CaitPolicies{k,1}(idx(subrange+1:end))=NaN;
                elseif strcmp(MF_CaitStrat_Order,'last')
                    IDX_CaitPolicies{k,1}(idx(1:end-subrange))=NaN;
                elseif strcmp(MF_CaitStrat_Order,'before')
                    IDX_CaitPolicies{k,1}(idx)=NaN; %Nan all labeled values since we are interested only in trials prior
                    if idx(1) > subrange
                        IDX_CaitPolicies{k,1}(idx(1)-subrange:idx(1)-1)=i; %Set the range prior to the current strat label
                    elseif idx(1) > 1
                        IDX_CaitPolicies{k,1}(1:idx(1)-1)=i; %Set the range prior to the current strat label
                    end
                elseif strcmp(MF_CaitStrat_Order,'after')
                    IDX_CaitPolicies{k,1}(idx)=NaN; %Nan all labeled values since we are interested only in trials after
                    if length(IDX_CaitPolicies{k,1})-(idx(end)+1) > subrange
                        IDX_CaitPolicies{k,1}(idx(end)+1:idx(end)+subrange)=i; %Set the range after to the current strat label
                    elseif idx(end)+1 < length(IDX_CaitPolicies{k,1})
                        IDX_CaitPolicies{k,1}(idx(end)+1:end)=i; %Set the range after to end of behavior to the current strat label
                    end
                end
            end

            idx = find(IDX_CaitPolicies{k,2}==i);
            if length(idx) > subrange
                if strcmp(MF_CaitStrat_Order,'first')
                    IDX_CaitPolicies{k,2}(idx(subrange+1:end))=NaN;
                elseif strcmp(MF_CaitStrat_Order,'last')
                    IDX_CaitPolicies{k,2}(idx(1:end-subrange))=NaN;
                elseif strcmp(MF_CaitStrat_Order,'before')
                    IDX_CaitPolicies{k,2}(idx)=NaN; %Nan all labeled values since we are interested only in trials prior
                    if idx(1) > subrange
                        IDX_CaitPolicies{k,2}(idx(1)-subrange:idx(1)-1)=i; %Set the range prior to the current strat label
                    elseif idx(1) > 1
                        IDX_CaitPolicies{k,2}(1:idx(1)-1)=i; %Set the range prior to the current strat label
                    end
                end
            end
        end
    end

    IDX_CaitCombined1{k,1} = nan(size(timepoints{k,2},1),1);
    IDX_CaitCombined1{k,2} = nan(size(timepoints{k,3},1),1);

    for tr=1:size(timepoints{k,2},1)
        if IDX_CaitClass{k,1}(tr)==2 || IDX_CaitClass{k,1}(tr)==3 %Divergent S+E
            IDX_CaitCombined1{k,1}(tr) = 2;
        elseif IDX_CaitPolicies{k,1}(tr)==1 %IA-Policy
            IDX_CaitCombined1{k,1}(tr) = 1;
        elseif IDX_CaitClass{k,1}(tr)==0 %Exploratory
            IDX_CaitCombined1{k,1}(tr) = 3;
        elseif IDX_CaitPolicies{k,1}(tr)==3 %RS-Policy
            IDX_CaitCombined1{k,1}(tr) = 4;
        end
    end

    for tr=1:size(timepoints{k,3},1)
        if IDX_CaitClass{k,2}(tr)==2 || IDX_CaitClass{k,2}(tr)==3 %Divergent S+E
            IDX_CaitCombined1{k,2}(tr) = 2;
        elseif IDX_CaitPolicies{k,2}(tr)==1 %IA-Policy
            IDX_CaitCombined1{k,2}(tr) = 1;
        elseif IDX_CaitPolicies{k,2}(tr)==3 %RS-Policy
            IDX_CaitCombined1{k,2}(tr) = 4;
        elseif IDX_CaitClass{k,2}(tr)==0 %Exploratory
            IDX_CaitCombined1{k,2}(tr) = 3;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IDX_CaitCombined2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1 = IA-Policy, 2 = Divergence, 3 = Exploratory, 4 = Acquisition, 5 = RS-Policy

    IDX_CaitCombined2{k,1} = nan(size(timepoints{k,2},1),1);
    IDX_CaitCombined2{k,2} = nan(size(timepoints{k,3},1),1);

    %%%%%IA Portion of Task%%%%%%
    for tr=1:size(timepoints{k,2},1)
        if IDX_CaitPolicies{k,1}(tr)==1 %IA-Policy
            IDX_CaitCombined2{k,1}(tr) = 1;
        elseif IDX_CaitClass{k,1}(tr)==0 %Exploratory
            IDX_CaitCombined2{k,1}(tr) = 3;
        elseif IDX_CaitPolicies{k,1}(tr)==3 %RS-Policy
            IDX_CaitCombined2{k,1}(tr) = 5;
        end
    end

    idx2 = find(IDX_CaitPolicies{k,1}==1); %Find all IA-Policy Stable Policies
    if ~isempty(idx2)
        idx2 = idx2([1;flip(diff(flip(idx2))~=-1)]==true); %Remove consecutive instances - e.g. keep only the 1st trial of a block of IA-Policy trials
        if idx2(end) > 1 %As long as there is a IA-Policy trial and it's not the 1st trial in RS...
            IDX_CaitCombined2{k,1}(idx2(end)-1)=4; %Set the t-1 of the 1st trial in the last IA-Policy block as the 'Acquisition' trial
        end
    else
        IDX_CaitCombined2{k,1}(1)=3; %Otherwise, animal acquired the IA-Policy on the first trial
    end

    %%%%%RS Portion%%%%%%
    for tr=1:size(timepoints{k,3},1)
        if IDX_CaitPolicies{k,2}(tr)==1 %IA-Policy
            IDX_CaitCombined2{k,2}(tr) = 1;
        elseif IDX_CaitPolicies{k,2}(tr)==3 %RS-Policy
            IDX_CaitCombined2{k,2}(tr) = 5;
            %elseif IDX_CaitClass{k,2}(tr)==0 %Exploratory
            %   IDX_CaitCombined2{k,2}(tr) = 3;
        end
    end

    clear idx1 idx2 idx3
    idx1 = find(IDX_CaitPolicies{k,2}==1); %Find all IA-Policy Stable Policies
    if ~isempty(idx1)%As long as there is an IA-Policy trial
        if idx1(end) < size(IDX_CaitPolicies{k,2},1)
            IDX_CaitCombined2{k,2}(idx1(end)+1)=2; %Set the t+1 of the last IA-Policy trial as the 'Divergent' trial
        end
    end

    idx2 = find(IDX_CaitPolicies{k,2}==3); %Find all RS-Policy Stable Policies
    idx2 = idx2(idx2>size(IDX_CaitPolicies{k,2},1)-10); %Remove Stable RS-Policies that were far away (>10 trials) from actual criterion

    if ~isempty(idx2)
        idx2 = idx2([1;flip(diff(flip(idx2))~=-1)]==true); %Remove consecutive instances - e.g. keep only the 1st trial of a block of RS-Policy trials
        if idx2(end) > 1 %As long as there is a RS-Policy trial and it's somehow not the 1st trial in RS...

            %IDX_CaitCombined2{k,2}(idx2(end)-1)=4; %Set the t-1 of the 1st trial in the last RS-Policy block as the 'Acquisition' trial

            %Optional: Make Acquisition the last error prior to the RS-Policy
            idx3 = find(cell2mat(timepoints{k,3}(:,6))==1); %Grab all errors
            idx3 = idx3(idx3 < idx2(end)); %Select all errors that precede the final RS-Policy block

            %if IDX_CaitPolicies{k,2}(idx3(end))~=1 %If excluding IA-Policy trials
            IDX_CaitCombined2{k,2}(idx3(end))=4; %Set the final error preceding the final RS-Policy block (within 10 trials of end of behavior) as acquisition
            %end
        end
    else
        idx2 = find(cell2mat(timepoints{k,3}(:,6))==1,1,'last');
        %IDX_CaitCombined2{k,2}(idx2)=4; %If the animal never technically read RS-Policy criterion, just set final error as acquisition
    end

    idx1 = find(IDX_CaitCombined2{k,2}==2);
    idx2 = find(IDX_CaitCombined2{k,2}==4);
    if ~isempty(idx1) && ~isempty(idx2) %TEST (commented out lines 372/373) - make exploratory trials between the IA-Divergence and RS-Acquisition
        if idx1(end) < idx2(end)
            IDX_CaitCombined2{k,2}(idx1(end)+1:idx2(end)-1) = 3;
        end
    end
end

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'CaitClass';
IDX_set{2,IDX_counter} = IDX_CaitClass;

if MF_CaitClass_subsetON == true
    stringAddon = sprintf(' - %s %d trials',MF_CaitClass_Order,MF_CaitClass_Range);
else
    stringAddon = "";
end
if MF_CaitClass_CmbnDiv == true
    IDX_set{3,IDX_counter} = {sprintf('Exploratory%s',stringAddon),sprintf('Stable%s',stringAddon),sprintf('Divergent%s',stringAddon)};
    IDX_set{4,IDX_counter} = [0:2];
else
    IDX_set{3,IDX_counter} = {sprintf('Exploratory%s',stringAddon),sprintf('Stable%s',stringAddon),sprintf('Div S-E%s',stringAddon),sprintf('Div S-S%s',stringAddon),sprintf('Lapse%s',stringAddon),sprintf('Div E-S%s',stringAddon)};
    IDX_set{4,IDX_counter} = [0:5];
end
IDX_set{5,IDX_counter} = 'Caitriona''s ML Classification method organizing trials into: Exploratory, Stable, or Divergent trials';

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'CaitClass: Stable Policies';
IDX_set{2,IDX_counter} = IDX_CaitPolicies;
IDX_set{3,IDX_counter} = {'Unstable','IA_Cue','!IA_Cue','RS_Cue','!RS_Cue','Left','Right'};
IDX_set{4,IDX_counter} = [0:6];
IDX_set{5,IDX_counter} = 'Further defines CaitClass'' stable policies into specific types';

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'CaitClass: IA-Policy + Div + Exp + RS-Policy';
IDX_set{2,IDX_counter} = IDX_CaitCombined1;
IDX_set{3,IDX_counter} = {'IA-Policy','Div (->S+E)','Exploratory','RS-Policy'};
IDX_set{4,IDX_counter} = [1:4];
IDX_set{5,IDX_counter} = 'Limits/Combines trial types to: IA-Policy, Divergent, Exploratory, and RS-Policy';

IDX_counter=IDX_counter+1;
IDX_set{1,IDX_counter} = 'CaitClass: IA-Policy + IA-Div + Exp + RS-Acquisition + RS-Policy';
IDX_set{2,IDX_counter} = IDX_CaitCombined2;
IDX_set{3,IDX_counter} = {'IA-Policy','IA-Divergent','Exploratory','RS-Acquisition','RS-Policy'};
IDX_set{4,IDX_counter} = [1:5];
IDX_set{5,IDX_counter} = 'Further limites trials to: IA-Policy, IA-Divergent, Exploratory, RS-Acquisition, and RS-Policy';



clear MF_fs MF_range MF_winRange MF_radRange MF_degRange MF_taskstg MF_taskwin MF_behWindow
clear rng_prestrk rng_poststrk rng_strkstart rng_peristrk rng_peristrat rng_IAEnd rng_IAStart
clear MF_CaitClass_CmbnDiv MF_CaitClass_subsetON MF_CaitClass_Range MF_CaitClass_Order MF_CaitStrat_subsetON MF_CaitStrat_Range MF_CaitStrat_Order

end