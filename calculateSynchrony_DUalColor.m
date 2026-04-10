%{
XCorr Version 3.PS.Dual - AJP - Sohal Lab - January 2021
%}

cleaningON = true; %Perform signal cleaning or not

frequency = 40; % frequency to analyze (in Hz)
bandwidth = (50-30) / frequency; % filter around frequency using this
%To calculate use: [Freq Range] / [Center Freq]
%eg: [50hz-30hz] / [40hz] = 0.5 bandwidth

window = round(data.x.fs); %fs is sampling frequencies. window is the number of samples in each second (window is 1sec)
dt = 1/data.x.fs; %timestep

nshuffles = 1; %number of shuffles

clear x z r shf95p shf5p shf95p2 shf5p2 cor1 cor2 corTD Yr coeff;
%convert datatype to double with sev function
x(1,:) = double(data.x.data(1,:)); %Left ACE-mNeon
x(2,:) = double(data.x.data(2,:)); %Left Varnam
x(3,:) = double(data.x.data(3,:)); %Right ACE-mNeon
x(4,:) = double(data.x.data(4,:)); %Right Varnam
x(5,:) = double(data.x.data(5,:)); %Left cyOFP
x(6,:) = double(data.x.data(6,:)); %Right cyOFP

N = length(x(1,:)); %length of data points there are

if cleaningON == true

    %First filter does a broadband filter
    f = 40; bw = 1.92;%f = 40; bw = 1.95;
    L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter **From WIKI: FIR filters express each output sample as a weighted sum of the last N input samples, where N is the order of the filter.
    B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter

    for i=1:6,
        r(i,:) = filtfilt(B, 1, x(i,:)); %r is filtered version of x (runs forward and backward using filtfilt -- if it goes in one direction, there will be a phase delay
    end

    %Next perform robust linear regression to find the best fit of tdTomato
    %onto mNeon and then substract
    warning('off','stats:statrobustfit:IterationLimit');
    numPeriods = 10;
    regressWin = round(window/(f*(1-0.5*bw))*numPeriods);  %Window to apply to signals for regression. Window is dependent on lowest frequency, for the broadband freq (1-79Hz) with optimal # of periods = 10.
    r = tempocleanNorm_2C(r,regressWin); %tempoclean performs robust linear regression on each hemisphere to clean mNeons
    warning('on','stats:statrobustfit:IterationLimit');

    %Second filter around frequency of interest
    f = frequency;
    bw = bandwidth;
    L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter
    B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter

    for i=1:6,
        z(i,:) = filtfilt(B, 1, r(i,:)); %z is filtered version of x (runs forward and backward using filtfilt -- if it goes in one direction, there will be a phase delay
    end

    %Next perform a final signal cleaning, minimally done here but sometimes
    %intrahemisphere noise persists and this removes it
    warning('off','stats:statrobustfit:IterationLimit');
    numPeriods = 10;
    regressWin = round(window/(f*(1-0.5*bw))*numPeriods);
    z = tempocleanNorm_2C(z,regressWin);
    warning('on','stats:statrobustfit:IterationLimit');

else %IF no signal cleaning - just apply a single filter around frequency of interest
    f = frequency;
    bw = bandwidth;
    L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter
    B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter


    for i=1:6,
        z(i,:) = filtfilt(B, 1, x(i,:)); %z is filtered version of x (runs forward and backward using filtfilt -- if it goes in one direction, there will be a phase delay
    end
end

  
%Finally run through L and R mNeons and compare via cross-correlation.
%In addition compare to top 95% shuffled data and bottom 5% shuffled
%data

n = 0;%going through vector from first timepoint to end. take data in 1sec chunks. n is chunk number
xCorrWin = floor(window/4); %Length of time (in samples) to perform cross correlations on
xCorrOvlp = 0;
xCorrStep = floor(window/4);

totalSS = 360; %Total amount of shift (goes half this left shifted (-) and half right shifted (+))
numSS = 16; %Total steps to get to shift. Will split half + and - around 0 degrees
xCorrShft = round(((window/f)/(360/totalSS))/numSS);

%Note: Channel 3 (Right mNeon) and 4 (Right tdTomato) gets shifted left (-) or right (+)
R = ceil((N-xCorrWin-(2*xCorrOvlp))*rand(1,nshuffles)); %Use this to get a random starting point for 100 shuffles (less shuffled)

for in=xCorrWin:xCorrStep:(N-xCorrWin-xCorrOvlp-(xCorrShft*numSS))
    
    n=n+1;
    
    %mNeon Signals
    cor1(1+numSS,n,:) = xcorr(z(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(3,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),0,'normalized'); %Perform cross correlation between both adjusted and filtered mNeons
    
    for j=1:numSS
        cor1(numSS-j+1,n,:) = xcorr(z(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(3,in-xCorrOvlp-(j*xCorrShft):in+xCorrWin+xCorrOvlp-(j*xCorrShft)),0,'normalized'); %Perform cross correlation with a phase shift
    end
    
    for j=1:numSS
        cor1(j+1+numSS,n,:) = xcorr(z(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(3,in-xCorrOvlp+(j*xCorrShft):in+xCorrWin+xCorrOvlp+(j*xCorrShft)),0,'normalized'); %Perform cross correlation with a phase shift
    end
    
    %Ruby Signals
    cor2(1+numSS,n,:) = xcorr(z(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(4,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),0,'normalized'); %Perform cross correlation between both adjusted and filtered mNeons
    
    for j=1:numSS
        cor2(numSS-j+1,n,:) = xcorr(z(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(4,in-xCorrOvlp-(j*xCorrShft):in+xCorrWin+xCorrOvlp-(j*xCorrShft)),0,'normalized'); %Perform cross correlation with a phase shift
    end
    
    for j=1:numSS
        cor2(j+1+numSS,n,:) = xcorr(z(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(4,in-xCorrOvlp+(j*xCorrShft):in+xCorrWin+xCorrOvlp+(j*xCorrShft)),0,'normalized'); %Perform cross correlation with a phase shift
    end

    corTD(1,n,:) = xcorr(z(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(6,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),0,'normalized'); %Perform cross correlation between left mNeon and right cyOFP - for investigating level of cyOFP contamination still present
    corTD(2,n,:) = xcorr(z(2,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(6,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),0,'normalized'); %Also do left varnam and right cyOFP
  
end


clear numShftSteps window tmp regressWin R R10 nshuffles numPeriods N n in i L B bw f data dt




            