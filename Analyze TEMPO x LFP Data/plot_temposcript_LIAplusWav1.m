%{
Modified version of plot_temposcript that cycles through a folder of data
structs exported from Adam's Open_ephys_loader_test1 app


%}


% Get a list of all files and folders in this folder.
files = dir('C:\Users\aphensy\Box\Postdoctoral Work\Gamma Oscillation Project\TEMPO Data\LFP + PVI - Carl and Kathleen\carl_lfp+tempo\Day 1 Tanks');
files(ismember( {files.name}, {'.', '..'})) = [];  %remove . and ..
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);


%Perform temposcript on each file - NOTE, only have the data structs in
%selected folder
clear sigMatrix
for k = 1 : length(files)
    
    data = SEV2mat(strcat(subFolders(k).folder,'\',subFolders(k).name));
    data = combineLIAWav1(data);
    
    temposcript_openephys
    
    sigMatrix{k,1} = files(k).name;
    sigMatrix{k,2} = strcat(num2str(frequency*(1-0.5*bandwidth)),'hz- ', num2str(frequency*(1+0.5*bandwidth)),'hz');
    
    %For standard xCorr

    sigMatrix{k,3} = cor1;
    sigMatrix{k,4} = lfpz1;
    sigMatrix{k,5} = txl1;
    sigMatrix{k,6} = pwr1;
    sigMatrix{k,7} = shf1;
    sigMatrix{k,8} = plv1;
    sigMatrix{k,9} = aec1;
    sigMatrix{k,10} = lag1;
    sigMatrix{k,11} = lfp1;
    sigMatrix{k,12}  = wpl1;
    sigMatrix{k,13} = shf2;
    sigMatrix{k,14} = gba1;
    
end
%% Set the k animal set and behavior variables initially

animalSet = [1:4];
sel = 1;
selections = {"Baseline","Initial-Association","Rule-Shift","X Cropped Minutes Before End"};
%transitions = round([20 35201 50751 42880; 20 25841 92460 92460; 20 42561 86651 85820]/20)*4; %BL->IA and IA->RS; animal 2 missing IA-RS transition and failed;  FPS = 20
transitions = round([20 470 1860 5446; 20 325 2141 4945; 20 404 3769 6233; 20 421 2230 4644])*10; %BL->IA and IA->RS; animal 2 missing IA-RS transition and failed;  FPS = 20
cropbehavior = false;

% Create a Z-Scored sigMatrix across all metrics
sigMatrix_Z = sigMatrix(animalSet,:);
for k=animalSet
    for i=3:size(sigMatrix,2)-1
        if cropbehavior
            sigMatrix_Z{k,i} = zscore(sigMatrix{k,i}(:,transitions(k,sel):end,:),0,2);
        else
            sigMatrix_Z{k,i} = zscore(sigMatrix{k,i}(:,1:end,:),0,2);
        end
    end
end

% Generate sigSHFComp based on the shf1 - Window x Window Threshold
clear sigSHFComp
AECThresh = prctile(sigMatrix_Z{k,9},90);
for k=animalSet
    n=0;
    if cropbehavior
        for i=transitions(k,sel):length(sigMatrix{k,3})
            n = n+1;
            sigSHFComp{k,1}(n) = sigMatrix{k,3}(1,i) > prctile(sigMatrix{k,7}(:,i), 95);
            sigSHFComp{k,2}(n) = sigMatrix{k,3}(1,i) < prctile(sigMatrix{k,7}(:,i), 5);
            sigSHFComp{k,3}(n) = sigMatrix{k,4}(1,i) > prctile(sigMatrix{k,13}(:,i), 95);
            sigSHFComp{k,4}(n) = (sigMatrix_Z{k,9}(1,n)) > AECThresh;
            sigSHFComp{k,5}(1,n) = logical(sigMatrix{k,14}(1,i));
            sigSHFComp{k,5}(2,n) = logical(sigMatrix{k,14}(2,i));
        end
    else
        for i=1:length(sigMatrix{k,3})
            n = n+1;
            sigSHFComp{k,1}(n) = sigMatrix{k,3}(1,i) > prctile(sigMatrix{k,7}(:,i), 95);
            sigSHFComp{k,2}(n) = sigMatrix{k,3}(1,i) < prctile(sigMatrix{k,7}(:,i), 5);
            sigSHFComp{k,3}(n) = sigMatrix{k,4}(1,i) > prctile(sigMatrix{k,13}(:,i), 95);
            sigSHFComp{k,4}(n) = (sigMatrix_Z{k,9}(1,i)) > AECThresh;
            sigSHFComp{k,5}(1,n) = logical(sigMatrix{k,14}(1,i));
            sigSHFComp{k,5}(2,n) = logical(sigMatrix{k,14}(2,i));
        end
    end
end

sigMatrix_Z_bkp = sigMatrix_Z;
sigSHFComp_bkp = sigSHFComp;
%% Create the z-score sigMatrix only
sigMatrix_Z = sigMatrix;
for k=1:size(sigMatrix,1)
    for i=3:size(sigMatrix,2)
        sigMatrix_Z{k,i} = zscore(sigMatrix{k,i}(:,1:end,:),0,2);
    end
end

%% Generate sigSHFComp based on the shf1 - GLOBAL THRESHOLD

for k=1:size(sigMatrix,1)
    %First take the maximum value across all windows per shuffle
    for s=1:size(sigMatrix{k,7},1)
        maxNull(s) = max(sigMatrix{k,7}(s,:));
    end
    %Then set the global threshold value based on percentile
    thr = prctile(sigMatrix{k,7}, 95,'all');
    %thr = prctile(maxNull, 95);
    
    n=0;
    for i=1:length(sigMatrix{k,3})
        n = n+1;
        sigSHFComp{k,1}(n) = sigMatrix{k,3}(1,i) > thr;
    end
end

%%
data = combineLIAWav1(data);

%% SUPPORT FUNCTIONS

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

data_out = struct('lfp_data',lfpData,'tempo_data',tempoData,'fs',fs);

end


%EOF


