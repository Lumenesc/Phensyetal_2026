function [brstCorr,brstCorr2,bursts] = calcBurstCorrelations(data,tempoFreq,range1,lfpFreq,range2,xCorrWin,maxLagtxl,stdDevThrsh,minBurstDur)

% %%%%%%%%%%%%%%%%calcBurstCorrelations %%%%%%%%%%%%%%%%
%               
%               %%%%%%%% OUTPUTS %%%%%%%%
% brstCorr: Correlation between the lfp and tempo signals within one
% recording site (concatenating site 1 and site 2) for each window which
% contained a burst
% 
% bursts: the onset (clmn1) and end (clmn2) of each identified burst, in
% seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



window = round(data.fs); %fs is sampling frequencies. window is the number of samples in each second (window is 1sec)
dt = 1/data.fs; %timestep
tempoBW = range1/tempoFreq;
lfpBW = range2/lfpFreq;

%Convert values from seconds to samples
xCorrWin = round(xCorrWin*window);
maxLagtxl = round(maxLagtxl*window);
minBurstDur = round(minBurstDur*window);

%convert datatype to double with sev function
x(1,:) = double(data.tempo_data(1,:)); %Signal Site 1
x(2,:) = double(data.tempo_data(2,:)); %Reference Site 1
x(3,:) = double(data.tempo_data(3,:)); %Signal Site 2
x(4,:) = double(data.tempo_data(4,:)); %Reference Site 2

lfp(1,:) = double(data.lfp_data(1,:)); %Signal Site 1
lfp(2,:) = double(data.lfp_data(2,:)); %Signal Site 1

%First filter does a broadband filter
f = 40; bw = 1.95;
L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter **From WIKI: FIR filters express each output sample as a weighted sum of the last N input samples, where N is the order of the filter.
B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter
B = B'; 
for i=1:4,
    r(i,:) = filtfilt(B, 1, x(i,:)); %r is filtered version of x (runs forward and backward using filtfilt -- if it goes in one direction, there will be a phase delay
end
%Next perform robust linear regression to find the best fit of tdTomato 
%onto mNeon and then substract
fprintf('First broadband cleaning...\n');
numPeriods = 10;
regressWin = round(window/(f*(1-0.5*bw))*numPeriods);  %Window to apply to signals for regression. Window is dependent on lowest frequency, for the broadband freq (1-79Hz) with optimal # of periods = 10.
r = tempocleanNorm(r,regressWin); %tempoclean performs robust linear regression on each hemisphere to clean mNeons

%Second filter around frequency of interest
f = tempoFreq;
bw = tempoBW;
L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter
B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter 
for i=1:4,
    z(i,:) = filtfilt(B, 1, r(i,:)); %z is filtered version of x (runs forward and backward using filtfilt -- if it goes in one direction, there will be a phase delay
end
    
%Next perform a final signal cleaning, minimally done here but sometimes 
%intrahemisphere noise persists and this removes it
fprintf('Second narrowband cleaning (%dhz - %dhz)...\n',tempoFreq*(1-0.5*tempoBW),tempoFreq*(1+0.5*tempoBW));
regressWin = round(window/(f*(1-0.5*bw))*numPeriods); 
z = tempocleanNorm(z,regressWin);

%Filter LFP at 30-80hz, apply a notch filter from 59:61hz first
f = lfpFreq;
bw = lfpBW;
L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter
B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter
for i=1:2
    if f*(1-0.5*bw) < 60 && f*(1+0.5*bw) > 60 %60hz noise is within filter range
        lfpZ(i,:) = bandstop(lfp(i,:),[59 61],window); %Perform notch filter from 59 to 61hz
        lfpZ(i,:) = filtfilt(B, 1, lfpZ(i,:)); %Also filter the LFP signals
    else
        lfpZ(i,:) = filtfilt(B, 1, lfp(i,:)); %Also filter the LFP signals
    end
end 


%Identify bursts for both Left and Right hemisphere - then concatenate lists
%envZ = zscore(abs(hilbert(lfpZ(1,:))));
envZ = calcZScoreEnv(lfpZ(1,:));
isHigh = envZ > stdDevThrsh;
isHigh(end-xCorrWin-1:end) = false; %Remove bursts that occur during the end of recording
burstsL = findBurstIndices(isHigh,minBurstDur);
%envZ = zscore(abs(hilbert(lfpZ(2,:))));
envZ = calcZScoreEnv(lfpZ(2,:));
isHigh = envZ > stdDevThrsh;
isHigh(end-xCorrWin-1:end) = false; %Remove bursts that occur during the end of recording
burstsR = findBurstIndices(isHigh,minBurstDur);

bursts = [burstsL;burstsR];
%Also identify random false bursts periods to use as a control
burstsRnd = double.empty(length(bursts),0);
burstsRnd(:,1)=randperm(length(lfpZ)-xCorrWin,length(bursts))';
burstsRnd(:,2)=randperm(length(lfpZ)-xCorrWin,length(bursts))';
%% Perform cross-correlation between LFP and TEMPO Signal during Bursts
brstCorr = cell(2,3);
brstCorr{1,1} = 'Bursts'; brstCorr{1,2} = 'False Bursts'; brstCorr{1,3} = 'Shuffled';
brstCorr(2,1:3) = {double.empty((maxLagtxl*2)+1,0)};
brstCorr2 = brstCorr;
for brst = 1:length(bursts) %%%%%%% Measure cross-correlation along length of xCorrWin starting at burst onset
    %Burst Onset
    brstCorr{2,1}(:,brst) = xcorr(normalize(z(1,bursts(brst,1):bursts(brst,1)+xCorrWin),"range",[-1,1]),normalize(lfpZ(1,bursts(brst,1):bursts(brst,1)+xCorrWin),"range",[-1,1]),maxLagtxl,'normalized')'; %Perform cross correlation between left mNeon and left LFP
    brstCorr{2,2}(:,brst) = xcorr(normalize(z(1,burstsRnd(brst,1):burstsRnd(brst,1)+xCorrWin),"range",[-1,1]),normalize(lfpZ(1,burstsRnd(brst,1):burstsRnd(brst,1)+xCorrWin),"range",[-1,1]),maxLagtxl,'normalized')'; %Perform cross correlation between left mNeon and left LFP
    brstCorr{2,3}(:,brst) = xcorr(normalize(z(1,burstsRnd(brst,1):burstsRnd(brst,1)+xCorrWin),"range",[-1,1]),normalize(lfpZ(1,burstsRnd(brst,2):burstsRnd(brst,2)+xCorrWin),"range",[-1,1]),maxLagtxl,'normalized')'; %Perform cross correlation between left mNeon and left LFP
end

for brst = 1:length(bursts) %%%%%%% Measure cross-correlation along length of the burst instead of a preset window
    %Burst Onset
    burstlength = bursts(brst,2)-bursts(brst,1);
    brstCorr2{2,1}(:,brst) = xcorr(normalize(z(1,bursts(brst,1):bursts(brst,2)),"range",[-1,1]),normalize(lfpZ(1,bursts(brst,1):bursts(brst,2)),"range",[-1,1]),maxLagtxl,'normalized')'; %Perform cross correlation between left mNeon and left LFP
    brstCorr2{2,2}(:,brst) = xcorr(normalize(z(1,burstsRnd(brst,1):burstsRnd(brst,1)+burstlength),"range",[-1,1]),normalize(lfpZ(1,burstsRnd(brst,1):burstsRnd(brst,1)+burstlength),"range",[-1,1]),maxLagtxl,'normalized')'; %Perform cross correlation between left mNeon and left LFP
    brstCorr2{2,3}(:,brst) = xcorr(normalize(z(1,burstsRnd(brst,1):burstsRnd(brst,1)+burstlength),"range",[-1,1]),normalize(lfpZ(1,burstsRnd(brst,2):burstsRnd(brst,2)+burstlength),"range",[-1,1]),maxLagtxl,'normalized')'; %Perform cross correlation between left mNeon and left LFP
end


%Convert bursts to seconds (for downstream analysis using timepoints)
bursts = (bursts/window);


end

function envZ = calcZScoreEnv(sig)
env = abs(hilbert(sig));
envMed = median(env);
envMad = mad(env, 1);              % median absolute deviation
envZ   = (env - envMed) / (1.4826*envMad + eps);
end

function bursts = findBurstIndices(isHigh,minDur)

d = diff([0; isHigh(:); 0]);
starts = find(d == 1);
ends   = find(d == -1) - 1;

dur = ends - starts + 1;
keep = dur >= minDur;

burstStarts = starts(keep);
burstEnds   = ends(keep);

bursts = [burstStarts burstEnds];

end

