%% Step 1 - Build the PV Synchrony Gate - E.g. Which trials have high probability of PVI Post-Outcome Synchrony
clear stateSummary synchSummary stateNames IDX_TITLES stateTable synchTable

%Get the probabilities of PVI across all trial's pre-decision to set significant threshold
THRESH_PERCENT = 90;
cState = [3];
getNextTrial = false; pair_selection = NaN;
IDX_SELECTION = [2]; 
IDX_RANGES = {[2]};
threshold = prctile(getProbabilities(IDX_ALL_INFO,IDX_ALL_CAT,IDX_set,stateArrayFull,IDX_SELECTION,IDX_RANGES,cState,getNextTrial,pair_selection),THRESH_PERCENT);

%%
getNextTrial = false; pair_selection = 1;

%Sets the IDX_SELECTION and IDX_RANGES based on the selection
cState = [3];
IDX_SELECTION = [1 2]; 
IDX_RANGES = {[3],[2]};
probabilities_pv_bkp = getProbabilities(IDX_ALL_INFO,IDX_ALL_CAT,IDX_set,stateArrayFull,IDX_SELECTION,IDX_RANGES,cState,getNextTrial,pair_selection);

cState = [40 30 38 31];
IDX_SELECTION = [1 2]; 
IDX_RANGES = {[3],[1]};
probabilities_bkp = getProbabilities(IDX_ALL_INFO,IDX_ALL_CAT,IDX_set,stateArrayFull,IDX_SELECTION,IDX_RANGES,cState,getNextTrial,pair_selection);
        
% probabilities_pv = probabilities_pv_bkp(1:end-1);
% probabilities = probabilities_bkp(2:end,:);
probabilities_pv = probabilities_pv_bkp;
probabilities = probabilities_bkp;
%%

pvGate = probabilities_pv >= threshold;

% GENERATE SHUFFLED PVGATE MATRIX
nPerm = 1000;
N = numel(pvGate);

pvGate_shuff = zeros(N, nPerm);

for p = 1:nPerm
    pvGate_shuff(:,p) = pvGate(randperm(N));
end
pvGate_shuff = logical(pvGate_shuff);

% GENERATE SHUFFLED DISTRIBUTIONS AND CALCULATE SIGNIFICANCE THRESHOLDS
probMeans_shuff = nan(size(pvGate_shuff,2),4);
sigThresh = nan(2,4);
probabilitiesSplit = cell.empty(0,2);
for p = 1:size(pvGate_shuff,2)
    probMeans_shuff(p,1:4) = mean(probabilities(pvGate_shuff(:,p),:),1) ./ mean(probabilities(~pvGate_shuff(:,p),:),1);
end
probMeans = mean(probabilities(pvGate,:),1) ./ mean(probabilities(~pvGate,:),1);
sigThresh(1,:) = prctile(probMeans_shuff,95,1);
sigThresh(2,:) = prctile(probMeans_shuff,5,1);
probabilitiesSplit{1} = probabilities(pvGate,:);
probabilitiesSplit{2} = probabilities(~pvGate,:);

% CALCULATE STATS
fprintf('\nCalculating Statistics: Shuffled Label Permutation and Rank Sum\n')
stateNames = {'Contra-IP','Contra-AP','Ipsi-IP','Ipsi-AP'};
for i=1:4
if probMeans(i) > sigThresh(1,i)
    fprintf('%s is higher than 95%% of shuffled data. p(%s) = %f.\n',stateNames{i},stateNames{i},probMeans(i));
    figure
    histogram(probabilitiesSplit{1}(:,i),10,"Normalization","probability");
    hold on
    histogram(probabilitiesSplit{2}(:,i),10,"Normalization","probability");
elseif probMeans(i) < sigThresh(2,i)
    fprintf('%s is lower than 5%% of shuffled data. p(%s) = %f.\n',stateNames{i},stateNames{i},probMeans(i));
    figure
    histogram(probabilitiesSplit{1}(:,i),10,"Normalization","probability");
    hold on
    histogram(probabilitiesSplit{2}(:,i),10,"Normalization","probability");
end

[p,~] = ranksum(probabilitiesSplit{1}(:,i),probabilitiesSplit{2}(:,i));
if p < 0.3
    fprintf('%s has ranksum p-value of %f\n',stateNames{i},p)
end

end


%% SUPPORT FUNCTIONS


function probabilities = getProbabilities(IDX_ALL_INFO,IDX_ALL_CAT,IDX_set,stateArrayFull,IDX_SELECTION,IDX_RANGES,cState,getNextTrial,pair_selection)

numIter = max(IDX_ALL_INFO{2,3}{1,6});
IDX_SELECTION = [12 6 IDX_SELECTION]; %Adds the Trial ID as the first SELECTION paramater
IDX_RANGES = [{0} {0} IDX_RANGES];

clear statePost_Mean
for k=1:numIter
    IDX_RANGES{2} = k;
    fprintf('Performing aggregation %d/%d\n',k,numIter)

    pairIDX = zeros(max(max(IDX_ALL_CAT(12,IDX_ALL_CAT(6,:)==k))),1); %Generate an IDX for the pairs concatenated for full length of IA+RS
    pairIDX(end-length(IDX_set{2,16}{k,1})+1:end) = IDX_set{2,16}{k,1}; 
    numIter2 = max(max(IDX_ALL_CAT(12,IDX_ALL_CAT(6,:)==k)));
    IDX_ALL_INFO{2,3}{12} = [1:numIter2];
    fprintf('   Trials (%d): ',numIter2);
    
    for tr=1:numIter2

        pairIDX(tr) = 0;
        idx = find(pairIDX==pair_selection);
        
        if ~getNextTrial && ~ismember(1,IDX_SELECTION)
            IDX_RANGES{1} = tr;
            fprintf('PostOut%d ',tr);
        elseif getNextTrial && ~isempty(idx) && tr > numIter2-length(IDX_set{2,16}{k,1}) && IDX_RANGES{4} == 1%Get the next conflict or non-conflict trial's values - nxt pre-dec
            IDX_RANGES{1} = idx(1);
            fprintf('PRE-DEC FOR NEXT TRIAL. T: %d, T+x: %d\n',tr,idx(1));
        elseif getNextTrial && ~isempty(idx) && tr > numIter2-length(IDX_set{2,16}{k,1}) && IDX_RANGES{4} == 2%Get the next conflict or non-conflict trial's values - current post-out
            IDX_RANGES{1} = tr;
            fprintf('POST-OUT FOR TRIAL T. T: %d, T+x: %d\n',tr,idx(1));
        elseif ~getNextTrial  && tr < numIter2 && IDX_RANGES{4} == 1 %Get the current trial's values
            IDX_RANGES{1} = tr+1;
            fprintf('PreDec%d ',tr);
        elseif ~getNextTrial  && tr < numIter2 && IDX_RANGES{4} == 2 %Get the current trial's values
            IDX_RANGES{1} = tr;
            fprintf('PostOut%d ',tr);
        else
            IDX_RANGES{1} = 1;
        end

        T=combinations(IDX_RANGES{:});
        T=unique(T); T=standardizeMissing(T,99); %Remove duplicate rows and then make 99s NaN - note: unique does not work with NaN values, thus this step is required
        
        n=0;
        for i=1:size(T,1)
            n=n+1;
            values = rmmissing(T{i,:});
            combIDX = ones(1,size(IDX_ALL_CAT,2));
            for j=1:size(values,2)
                combIDX = combIDX & IDX_ALL_CAT(IDX_SELECTION(j),:)==IDX_ALL_INFO{2,3}{IDX_SELECTION(j)}(values(j));
            end
            statePost_Mean{k,1}{tr,n} = stateArrayFull(cState,combIDX)'; %Mean probability of the state existing in the behavioral windows
        end
    end
    fprintf('\n');
end

clear validTrials probabilities
probabilities = [];
n=0;
for k = 1:size(statePost_Mean,1)
    for con = 1:size(statePost_Mean{k,1},2)
        validTrials{k,1} = statePost_Mean{k,1}(~cellfun(@isempty,statePost_Mean{k,1}(:,con)),con);
        for i=1:length(validTrials{k,1})
            validTrials{k,1}{i,2} = mean(validTrials{k,1}{i,1}); %Post-Outcome Concatenation
        end
        probabilities = [probabilities;cell2mat(validTrials{k,1}(1:end-1,2))];
    end
end

end