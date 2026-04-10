function [] = findoscmodes(x, phaseoffsets, period, fs, win, inc, range,sigNames)

%%%%Aarron TEMPO assignments%%%%%

period = 1/40; %40Hz
fs = data.x.fs;
win = 0.250; %250ms
inc = 0.250; %250ms, no overlap
range = 1; %Examining at just a single cycle max

sigNames = {'PVL','MDL','PVR','MDR','RFL','RFR'};
if isempty(sigNames)
    for i=1:size(x,1)
        sigNames{i} = sprintf('Ch%d',i);
    end
end

% identify modes of oscillatory activity

% x = the inputs signals, assumed to be filtered in the band of interest

% phaseoffsets = the offset for each signal relative to the first (which is
% assumed to serve as the reference)

% period = period (in seconds) corresponding to frequency being analyzed

% fs = sampling frequency (per second)

% win = window size (in seconds) to use for calculating cross-correlogram

% inc = increment (in seconds) between successive windows

% range = max number of cycles to use for lag/lead, e.g., range = 1 means
% use one full cycle before / after

%%%%Aarron NOTE%%%%: Ch1 is Left-Hemisphere PVI

maxlag = round(range * period * fs); % number of indices of shift over
                                     % which to calculate cross-correlation

winsize = round(win*fs);
incsize = round(inc*fs);



% generate an array of synchrony values for each signal references to the
% first signal. Assume two possible phase differences, the preferred phase
% (specified by phaseoffsets) and another phase 180 degrees from that one

x = [masterx;zeros(2,1) masterxREF];
N = size(x);
nsig = N(1); % number of signals
npts = N(2); % length of signals
numIter = floor(npts/incsize)-ceil(winsize/incsize); %Number of iterations that will be needed to process the data

PVL=1; MDL=2; PVR=3; MDR=4; REFL=5; REFR=6;
sigCompSet = [[PVL,MDL,0,false,false,false];[PVL,PVR,0,false,false,false];[PVL,MDR,0,false,false,false];[PVR,MDL,0,false,false,false];[MDL,MDR,0,false,false,false];[PVR,MDR,0,false,false,false]]; %[Ch1,Ch2(will be shifted),phase-shift (in radians),includeLOWSynchstates,applyThreshold,mandatory state (NaN non-significant timepoints)] State table 10
shuffleON = false; %Will perform 100 shuffled data - NOTE will take a long time
thresholdON = false;
nshuffles = 100; perchigh = 90; perclow = 10; phaseThresh=85;
rndSeed = 1337;
plotpolarON = false;

smoothState = false;
smoothWindow = 4;

if isempty(gcp('nocreate'))
    parpool("threads");
end

close all

clear syncIn syncOut state phasesRawAll phasesRawAllThresh
if shuffleON == true
    clear thresh_syncIn thresh_syncOut
end

state = NaN(size(sigCompSet,1),numIter);

for s=1:size(sigCompSet,1)
    i = sigCompSet(s,1);
    j = sigCompSet(s,2);

    if shuffleON == true
        clear temp
        rng(rndSeed); %Provides a seed for reproducibility
        fprintf('Compare %s vs %s Shuffle: ',sigNames{i},sigNames{j});
        for shf=1:nshuffles
            fprintf('%d ',shf);
            n=0;
            C = [];
            xSHF = circshift(x(j,:),randi([1 npts],1,1));
            parfor n=1:numIter
                ind = (n*incsize):(n*incsize)+winsize; %indices from current point to length of win
                C(n,:) = xcorr(x(i,ind), xSHF(1,ind), maxlag, 'normalized'); %Perform x-correlation with phase-shifts (up to the maxlag)
            end

            [maxcorr, phasediff] = max(C'); %Find the maximum correlation - value (maxcorr) and index (phasediff) for each timepoint calculated above. Note: phasediff will be a positive integer from 1:2*maxlag+1 - however 1:maxlag will be in 1 direction and maxlag:2*maxlag+1 will be in the opposite
            maxcorr = maxcorr';
            phasediff = phasediff';
            phases = pi*(phasediff-(maxlag+1))/maxlag; % convert phasedifferences from indices to radians
            phases = mod(phases - sigCompSet(s,3), 2*pi); %mod(phase,2*pi) wraps the phase to 2*pi
            temp{1}(:,shf) = (maxcorr .* (phases <= pi/4 | phases >= 7*pi/4)); %Grabs the r-values from the xcorr between signals that are 'in-phase' (-45:45deg) of the offsetted phase
            temp{2}(:,shf) = (maxcorr .* (phases >= 3*pi/4 & phases <= 5*pi/4)); %Grabs the r-values from the xcorr between signals that are 'out-of-phase' (135:225deg) of the offsetted phase
        end
        fprintf('100%% - Done\n');

        thresh_syncIn{s,1} = nonzeros(temp{1});
        thresh_syncOut{s,1} = nonzeros(temp{2});


    end
    clear xSHF temp;

    %Process real data
    fprintf('Compare %s vs %s ',sigNames{i},sigNames{j});
    n=0;
    C = []; CSHF = [];  
    parfor n=1:numIter
        ind = (n*incsize):(n*incsize)+winsize; %indices from current point to length of win
        C(n,:) = xcorr(x(i,ind), x(j,ind), maxlag, 'normalized'); %Perform x-correlation with phase-shifts (up to the maxlag)
    end
    [maxcorr, phasediff] = max(C'); %Find the maximum correlation - value (maxcorr) and index (phasediff) for each timepoint calculated above. Note: phasediff will be a positive integer from 1:2*maxlag+1 - however 1:maxlag will be in 1 direction and maxlag:2*maxlag+1 will be in the opposite
    maxcorr = maxcorr';
    phasediff = phasediff';
    phasesRaw = pi*(phasediff-(maxlag+1))/maxlag; % convert phasedifferences from indices to radians
    phases = mod(phasesRaw - sigCompSet(s,3), 2*pi); %mod(phase,2*pi) wraps the phase to 2*pi
    syncIn(s,:) = maxcorr .* (phases <= pi/4 | phases >= 7*pi/4); %Grabs the r-values from the xcorr between signals that are 'in-phase' (-45:45deg) of the offsetted phase
    syncOut(s,:) = maxcorr .* (phases >= 3*pi/4 & phases <= 5*pi/4); %Grabs the r-values from the xcorr between signals that are 'out-of-phase' (135:225deg) of the offsetted phase
    if thresholdON == true
        clear threshold
        threshold(1,1) = prctile(thresh_syncIn{s,1},perchigh);
        threshold(1,2) = prctile(thresh_syncIn{s,1},perclow);
        threshold(2,1) = prctile(thresh_syncOut{s,1},perchigh);
        threshold(2,2) = prctile(thresh_syncOut{s,1},perclow);
        
        if sigCompSet(s,4) == true %Low synchrony states set to ON
            state(s,:) = 2*(maxcorr<threshold(1,2));
            if sigCompSet(s,5) == true %Threshold values set to ON
                state(s,:) = state(s,:)' + ((phases <= pi/4 | phases >= 7*pi/4) & (syncIn(s,:)>threshold(1,1))') + -((phases >= 3*pi/4 & phases <= 5*pi/4) & (syncOut(s,:)>threshold(2,1))');
            else
                state(s,:) = state(s,:)' + ((phases <= pi/4 | phases >= 7*pi/4) & (syncIn(s,:)>threshold(1,2))') + -((phases >= 3*pi/4 & phases <= 5*pi/4) & (syncOut(s,:)>threshold(2,2))');
            end
        else
            if sigCompSet(s,5) == true %Threshold values set to ON
                state(s,:) = ((phases <= pi/4 | phases >= 7*pi/4) & (syncIn(s,:)>threshold(1,1))') + -((phases >= 3*pi/4 & phases <= 5*pi/4) & (syncOut(s,:)>threshold(2,1))');
            else
                state(s,:) = (phases <= pi/4 | phases >= 7*pi/4) + -(phases >= 3*pi/4 & phases <= 5*pi/4);
            end
        end

    else
        state(s,:) = (phases <= pi/4 | phases >= 7*pi/4) + -(phases >= 3*pi/4 & phases <= 5*pi/4);
    end
    fprintf('100%% - Done\n');

    phasesRawThresh = phasesRaw; phasesRawThresh(maxcorr<prctile(thresh_syncIn{s,1},phaseThresh))=NaN; %Apply the treshold to phasesRaw to NaN any non-significant timepoints
    
    phasesRawAll(s,:) = phasesRaw; %Combine all phasesRaw into one matrix
    phasesRawAllThresh(s,:) = phasesRawThresh; %Combine all phasesRawThresh into one matrix

    if plotpolarON==true
        figure
        polarhistogram(phasesRaw(maxcorr>threshold(1,1)),'BinEdges',deg2rad([-7.5:15:352.5]));
        title(sprintf('%s vs %s',sigNames{i},sigNames{j}))
    end
    
end


if smoothState == true
    temp = state(s,:);
    for i=ceil(smoothWindow/2):size(state,2)-ceil(smoothWindow/2)
        if state(s,i)~=0
            for j=1:round(smoothWindow/2)
                if state(s,i-j)==0
                    temp(1,i-j) = state(s,i);
                end
                if state(s,i+j)==0
                    temp(1,i+j) = state(s,i);
                end
            end
        end
    end
    state(s,:) = temp;
end


%NaN points in which ANY thresholded states are non-significant
for s=1:size(sigCompSet,1)
    if sigCompSet(s,6) == true
        for i=1:size(state,2)
            if state(s,i)==0
                state(:,i)=NaN;
            end
        end
    end
end



binRange = 15; %In Degrees
if plotpolarON==true
    figure
    polarhistogram([phasesRawAllThresh(3,:) phasesRawAllThresh(4,:)],'BinEdges',deg2rad([-binRange/2:binRange:360-(binRange/2)]));
    title('All xH-PVxMD');
    figure
    polarhistogram([phasesRawAllThresh(1,:) phasesRawAllThresh(6,:)],'BinEdges',deg2rad([-binRange/2:binRange:360-(binRange/2)]));
    title('All wH-PVxMD');
    distFig
end

clear synclabels statelabels
m=1; n=0;
for s=1:size(sigCompSet,1)
    i = sigCompSet(s,1);
    j = sigCompSet(s,2);

        synclabels{m,1} = sprintf('IP %s vs %s',sigNames{i},sigNames{j});
        synclabels{m+1,1} = sprintf('OOP %s vs %s',sigNames{i},sigNames{j});
        m=m+2;

        n=n+1;
        statelabels{n} = sprintf('%sx%s',sigNames{i},sigNames{j});
end

stateBKP = state;

clear period fs range maxlag winsize incsize N nsig npts PVL MDL PVR MDR shuffleON thresholdON nshuffles perchigh perclow rndSeed x m n i j


%% OR OPERATION + IDENTIFICATION OF SINGLE wH AS SAME OR OPPOSITE SIDE WHEN SINGLE xH-SYNCH
state = stateBKP;
for i=1:size(state,2)

    if isnan(state(1,i))
        state(7:16,i) = NaN;
    else

        %Perform wH-LH vs wH-RH
        if state(1,i) == state(6,i) %wH-LHvRH-Same Phase (Or both 0)
            state(7,i) = state(1,i);
        elseif state(1,i) == 0 %LH has no synch, but RH does
            state(8,i) = state(6,i);
        elseif state(6,i) == 0 %RH has no synch, but LH does
            state(8,i) = state(1,i);

        elseif state(1,i)*state(6,i) == -1 %wH-LHvRH- Opposite Phase
            state(9,i) = -3;
        end

        %Perform xH-LH vs xH-RH
        if state(3,i) == state(4,i) %xH-LHvRH-Same Phase
            state(10,i) = state(3,i);
        elseif state(3,i) == 0 %LH has no synch, but RH does
            state(11,i) = state(4,i);

            if state(1,i) == 0 && state(6,i) ~= 0 %wH: LH has no synch, but RH does (same PV x MD wH has synch)
                state(13,i) = state(6,i);
            elseif state(1,i) ~=0 && state(6,i) == 0 %wH: RH has no synch, but LH does (opposite PV x MD wH has synch)
                state(14,i) = state(1,i);
            end

        elseif state(4,i) == 0 %RH has no synch, but LH does
            state(11,i) = state(3,i);
            
            if state(6,i) == 0 && state(1,i) ~= 0 %wH: RH has no synch, but LH does (same PV x MD wH has synch)
                state(13,i) = state(1,i);
            elseif state(6,i) ~=0 && state(1,i) == 0 %wH: LH has no synch, but RH does (opposite PV x MD wH has synch)
                state(14,i) = state(6,i);
            end

        elseif state(3,i)*state(4,i) == -1 %xH-LHvRH- Opposite Phase
            state(12,i) = -3;
        end

        if state(2,i) == 0 % PVI is neither in Phase or anti-phase
            state(15,i) = 1;
        end


    end
end

state = state([2 7:end],:); %REMOVE MDxMD FROM THE STATE DIAGRAMS!!

clear statelabels
statelabels{1} = "PVLxPVR";
statelabels{2} = "wH-Both";
statelabels{3} = "wH-Single";
statelabels{4} = "wH-Opp";
statelabels{5} = "xH-Both";
statelabels{6} = "xH-Single";
statelabels{7} = "xH-Opp";
statelabels{8} = "wHwhenxH-SamePV";
statelabels{9} = "wHwhenxH-SameMD";
statelabels{10} = "PV has no phase";

%% SECOND OPTION - COMBINE ALL SYMETRICAL PHASES - EITHER IN PHASE OR ANTI PHASE IPSILATERAL OR CONTRALATERAL
state = stateBKP;
for i=1:size(state,2)

    if isnan(state(1,i))
        state(7:8,i) = NaN;
    else

        %Perform wH-LH vs wH-RH
        if (state(1,i) == state(6,i) && state(1,i) ~= 0) || (state(3,i) == state(4,i) &&  state(3,i) ~= 0) %wH-LHvRH-Same Phase or wH-LHvRH-Same Phase (not 0)
            state(7,i) = 1;
        else
            state(7,i) = 0;
        end        

        if state(2,i) == 0 % PVI is neither in Phase or anti-phase
            state(8,i) = 1;
        end


    end
end
state = state([2 7:end],:); %REMOVE MDxMD FROM THE STATE DIAGRAMS!!


clear statelabels
statelabels{1} = "PVLxPVR";
statelabels{2} = "PVIxMD-Symmetrical";
statelabels{3} = "PV has no phase";


%%

[C,~,~] = unique(state','rows'); %Get all unique combinations that exist in the dataset
if any(isnan(C))
    C(isnan(C(:,1)),:) = []; % remove all nans
    C(end+1,:) = NaN; % add the unique one.
end
%Find unique instances of C, ignoring nested subsets
% clear icSum CSorted
% icSum = accumarray(ic,1); %Summate the # of each unique event
% [~,I] = sort(icSum,'descend');
% CSorted = [C(I,:) icSum(I)];

%Find all occurances of C - including nested subsets (e.g. [0 1 0] is a subset of [1 1 0] and thus [1 1 0] is included in the count for [0 1 0]

clear CSum CSorted stateArrayFull
stateArrayFull = ones(size(C,1),size(state,2));
for cstate=1:size(C,1)
    for i=1:size(state,2)
        if all(isnan(state(:,i)))
            stateArrayFull(cstate,i) = 0;
        elseif all(C(cstate,:)==0) && ~all(state(:,i)==0)
            stateArrayFull(cstate,i) = 0;
        else
            for j=1:length(C(cstate,:))
                if C(cstate,j)~=0 && C(cstate,j)~=state(j,i)
                    stateArrayFull(cstate,i) = 0;
                end
            end
        end
    end
end
CSum = sum(stateArrayFull,2);
[~,I] = sort(CSum,'descend');
CSorted = [C(I,:) CSum(I)];
stateArrayFull = stateArrayFull(I,:);

clear CText
for i=1:size(CSorted,1)
    for j=1:size(CSorted,2)-1
        if CSorted(i,j) == 1
            CText{i,j} = sprintf('%s IP',statelabels{j});
        elseif CSorted(i,j) == -1
            CText{i,j} = sprintf('%s OOP',statelabels{j});
        elseif CSorted(i,j) == 2
            CText{i,j} = sprintf('%s Low',statelabels{j});
        elseif CSorted(i,j) == -3
            CText{i,j} = sprintf('%s',statelabels{j});
        else
            CText{i,j} = '';
        end
    end
    format long
    CText{i,j+1} = sprintf('Count: %d',CSorted(i,j+1));
    CText{i,j+2} = sprintf('%.2f%%',round(CSorted(i,j+1)/size(state,2)*100,2));
end

clear stateArray %Create a stateArray that lists each timepoint as a state # instead of a binary sequence
maxStates = size(CSorted,1);
for i=1:size(state,2)
    [~,idx] = ismember(state(:,i)',CSorted(:,1:end-1),"rows");
    if idx <= maxStates
        stateArray(1,i) = idx;
        
        if sum(abs(state(:,i))) == 0
            stateArray(2,i) = 0;
        else
            stateArray(2,i) = (sum(syncIn(:,i))+sum(syncOut(:,i))) / (nnz(syncIn(:,i)>0)+nnz(syncOut(:,i)>0));
        end

    else
        stateArray(1,i) = NaN;
        stateArray(2,i) = NaN;
    end
end

clear i j ia ic I idx n maxStates