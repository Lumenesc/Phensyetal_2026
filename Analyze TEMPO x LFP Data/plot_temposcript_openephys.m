%{
Modified version of plot_temposcript that cycles through a folder of data
structs exported from Adam's Open_ephys_loader_test1 app


%}


% Get a list of all files and folders in this folder.
files = dir('C:\Users\AJPhe\Box\Postdoctoral Work\Gamma Oscillation Project\TEMPO Data\LFP x PFC-MD\Aligned TEMPO x LFP Data');
files(ismember( {files.name}, {'.', '..'})) = [];  %remove . and ..


%Perform temposcript on each file - NOTE, only have the data structs in
%selected folder
clear sigMatrix
for k = 1 : length(files)
    
    load([files(k).folder '\' files(k).name]);
    fprintf('Performing temposcript on data file #%d/%d = %s...\n', k, length(files), files(k).name);
    
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

animalSet = [1:3];
sel = 3;
selections = {"Baseline","Initial-Association","Rule-Shift","X Cropped Minutes Before End"};
transitions = round([20 35201 50751 42880; 20 25841 92460 92460; 20 42561 86651 85820]/20)*4; %BL->IA and IA->RS; animal 2 missing IA-RS transition and failed;  FPS = 20
cropbehavior = false;

% Create a Z-Scored sigMatrix across all metrics
sigMatrix_Z = sigMatrix;
for k=1:size(sigMatrix,1)
    for i=3:size(sigMatrix,2)
        if cropbehavior
            sigMatrix_Z{k,i} = zscore(sigMatrix{k,i}(:,transitions(k,sel):end,:),0,2);
        else
            sigMatrix_Z{k,i} = zscore(sigMatrix{k,i}(:,1:end,:),0,2);
        end
    end
end

% Generate sigSHFComp based on the shf1 - Window x Window Threshold
clear sigSHFComp
for k=1:size(sigMatrix,1)
    n=0;
    for i=transitions(k,sel):length(sigMatrix{k,3})
        n = n+1;
        if cropbehavior
            sigSHFComp{k,1}(n) = sigMatrix{k,3}(1,i) > prctile(sigMatrix{k,7}(:,i), 95);
            sigSHFComp{k,2}(n) = sigMatrix{k,3}(1,i) < prctile(sigMatrix{k,7}(:,i), 5);
            sigSHFComp{k,3}(n) = sigMatrix{k,4}(1,i) > prctile(sigMatrix{k,13}(:,i), 95);
        else
            sigSHFComp{k,1}(n) = sigMatrix{k,3}(1,i) > prctile(sigMatrix{k,7}(:,i), 95);
            sigSHFComp{k,2}(n) = sigMatrix{k,3}(1,i) < prctile(sigMatrix{k,7}(:,i), 5);
            sigSHFComp{k,3}(n) = sigMatrix{k,4}(1,i) > prctile(sigMatrix{k,13}(:,i), 95);
        end
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


%EOF