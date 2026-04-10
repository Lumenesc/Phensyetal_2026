%{

XCorr Version 3 - FOR OpenEphysLFP - AJP - Sohal Lab - Nov 2020

Adapted to work with the data struct that is exported from Adam's
Open_ephys_loader_test1 app

Note that the data structure organization is different than TDT's struct

%}

nshuffles = 100;

frequency = 40; % center frequency to analyze (in Hz)
range = 20; %total range of frequency to analyze (eg center freq of 40Hz and range of 20hz =  30-50Hz);

bandwidth = range/frequency;

window = round(data.fs); %fs is sampling frequencies. window is the number of samples in each second (window is 1sec)
dt = 1/data.fs; %timestep

clear x z r shf1 shf95p shf5p cor1 Yr coeff lfp lfp1 lfpz1 txl1 tmp lfpZ pwr1 plv1 aec1 lag1 wpl1 burstsL burstsR gba1;
%convert datatype to double with sev function
x(1,:) = double(data.tempo_data(1,:)); %Signal Site 1
x(2,:) = double(data.tempo_data(2,:)); %Reference Site 1
x(3,:) = double(data.tempo_data(3,:)); %Signal Site 2
x(4,:) = double(data.tempo_data(4,:)); %Reference Site 2

lfp(1,:) = double(data.lfp_data(1,:)); %Signal Site 1
lfp(2,:) = double(data.lfp_data(2,:)); %Signal Site 1

N = min(length(x(1,:)),length(lfp(1,:))); %length of data points there are

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
f = frequency;
bw = bandwidth;
L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter
B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter
    
for i=1:4,
    z(i,:) = filtfilt(B, 1, r(i,:)); %z is filtered version of x (runs forward and backward using filtfilt -- if it goes in one direction, there will be a phase delay
end
    
%Next perform a final signal cleaning, minimally done here but sometimes 
%intrahemisphere noise persists and this removes it
fprintf('Second narrowband cleaning (%dhz - %dhz)...\n',frequency*(1-0.5*bandwidth),frequency*(1+0.5*bandwidth));
regressWin = round(window/(f*(1-0.5*bw))*numPeriods); 
z = tempocleanNorm(z,regressWin);

%Filter LFP at 30-80hz, apply a notch filter from 59:61hz first
% f = 55;
% bw = 50/55;
% L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter
% B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter
for i=1:2
    %lfpZ(i,:) = bandstop(lfp(i,:),[59 61],window);
    lfpZ(i,:) = filtfilt(B, 1, lfp(i,:)); %Also filter the LFP signals
end 
  
%Finally run through L and R mNeons and compare via cross-correlation.
%In addition compare to top 95% shuffled data and bottom 5% shuffled
%data

n = 0;%going through vector from first timepoint to end. take data in 1sec chunks. n is chunk number
xCorrWin = round(window/10); %Length of time (in samples) to perform cross correlations on
xCorrOvlp = 0;
maxLagtxl = round (window/80);

minBurstDur = round((window/40)*2);

%envZ = zscore(abs(hilbert(lfpZ(1,:))));
envZ = calcZScoreEnv(lfpZ(1,:));
isHigh = envZ > 2.0;
burstsL = findBurstIndices(isHigh,minBurstDur);
%envZ = zscore(abs(hilbert(lfpZ(2,:))));
envZ = calcZScoreEnv(lfpZ(2,:));
isHigh = envZ > 2.0;
burstsR = findBurstIndices(isHigh,minBurstDur);

fprintf('Perform xcorr across recording. Window: %d Overlap: %d...\n',xCorrWin,xCorrOvlp);
%% Perform xCorrs across recording
R = ceil((N-xCorrWin-(2*xCorrOvlp))*rand(1,nshuffles)); %Use this to get a random starting point for 100 shuffles (less shuffled)
for in=xCorrWin:xCorrWin:(N-xCorrWin-xCorrOvlp),
        
    n=n+1;
    cor1(1,n,:) = xcorr(z(3,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),0,'normalized'); %Perform cross correlation between both adjusted and filtered mNeons
    cor1(2,n,:) = xcorr(z(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(4,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),0,'normalized'); %Perform cross correlation between left mNeon and right tdTomato - for investigating level of tdTomato contamination still present

    txl1(1,n,:) = max(xcorr(normalize(z(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),"range",[-1,1]),normalize(lfpZ(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),"range",[-1,1]),maxLagtxl,'normalized')); %Perform cross correlation between left mNeon and left LFP
    txl1(2,n,:) = max(xcorr(normalize(z(3,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),"range",[-1,1]),normalize(lfpZ(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),"range",[-1,1]),maxLagtxl,'normalized')); %Perform cross correlation between right mNeon and right LFP

    [lfp1(1,n,:), lag] = xcorr(lfp(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),lfp(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),round(0.02*window),'normalized'); %Perform cross correlation between two unfiltered LFP channels, allow for a 20ms lead/lag
    [lfpz1(1,n,:), lag] = xcorr(lfpZ(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),lfpZ(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),round(0.02*window),'normalized'); %Perform cross correlation between two filtered LFP channels, allow for a 20ms lead/lag
    lag1(1,n,:) = (lag - (round(0.02*window)+1)) * 1000/window; %Phase-lag of LFP signals

    phiL = angle(hilbert(lfpZ(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp)));
    phiR = angle(hilbert(lfpZ(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp)));
    plv1(1,n,:) = abs(mean(exp(1i*(phiL - phiR)))); %Phase-Locking Value

    pwr1(1,n,:) = bandpower(lfp(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp), window, [30 50]);
    pwr1(2,n,:) = bandpower(lfp(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp), window, [30 50]);

    envL = abs(hilbert(lfpZ(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp)));
    envR = abs(hilbert(lfpZ(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp)));
    aec1(1,n,:) = corr(envL(:), envR(:));  %Amplitude-Envelope Coupling: Correlation of 'gamma bursts'

    axL = hilbert(lfpZ(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp));
    axR = hilbert(lfpZ(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp));

    wpl1(1,n,:) = wpli_from_analytic(axL, axR);

    %Is there a gamma burst within this window on either left or right
    %hemisphere
    if any(in<burstsL(:,1) & in+xCorrWin>burstsL(:,1)) || any(in<burstsL(:,2) & in+xCorrWin>burstsL(:,2))
        gba1(1,n,:) = 1;
    else
        gba1(1,n,:) = 0;
    end
    if any(in<burstsR(:,1) & in+xCorrWin>burstsR(:,1)) || any(in<burstsR(:,2) & in+xCorrWin>burstsR(:,2))
        gba1(2,n,:) = 1;
    else
        gba1(2,n,:) = 0;
    end


    %R = ceil((N-xCorrWin-(2*xCorrOvlp))*rand(1,nshuffles)); %Use this to get a new random position every iterration (more shuffled)
    clear tmp
    for i=1:nshuffles, %perform xcorr on shuffled right signal
    shf1(i,n,:) = xcorr(z(3,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(1,R(i):R(i)+xCorrWin+(2*xCorrOvlp)),0,'normalized'); %Left GEVI vs shuffled Right GEVI
    shf2(i,n,:) = xcorr(lfpZ(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),lfpZ(2,R(i):R(i)+xCorrWin+(2*xCorrOvlp)),0,'normalized'); %Left LFP vs shuffled Right LFP
    end
          
    % shf95p(1,n) = prctile(tmp,95);
    % shf5p(1,n) = prctile(tmp,5);
       
end


clear xCorrWin xCorrOvlp window tmp regressWin R R10 nshuffles numPeriods N n in i L B bw f data dt


function w = wpli_from_analytic(ax, ay)
% ax, ay: analytic signals (complex), same length vectors
% returns scalar wPLI for that window

imX = imag(ax .* conj(ay));      % Imaginary component of cross-spectrum (time-domain proxy)
den = sum(abs(imX));

if den == 0
    w = NaN;                     % or 0, but NaN is safer
else
    w = abs(sum(imX)) / den;     % wPLI
end
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

