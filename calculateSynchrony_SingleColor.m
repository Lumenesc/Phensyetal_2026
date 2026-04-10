

%{
XCorr Version 3.PS - AJP - Sohal Lab - January 2021
%}

cleaningON = true; %Perform signal cleaning or not
shuffleON = false;

frequency = 40; % frequency to analyze (in Hz)
bandwidth = (50-30) / frequency; % filter around frequency using this
%To calculate use: [Freq Range] / [Center Freq]
%eg: [50hz-30hz] / [40hz] = 0.5 bandwidth

sigCHs = [1 3]; %What channels are the two signal carriers (e.g. mNeons)
ctrlCHs = [2 4]; %What channels to compare as a control check (e.g. mNeon CH1 vs tdTomato CH4)

window = round(data.x.fs); %fs is sampling frequencies. window is the number of samples in each second (window is 1sec)
dt = 1/data.x.fs; %timestep

nshuffles = 100; %number of shuffles

clear x z r shf95p shf5p cor1 corTD Yr coeff;
x(1,:) = double(data.x.data(1,:)); %convert datatype to double with sev function
x(2,:) = double(data.x.data(2,:));
x(3,:) = double(data.x.data(3,:));
x(4,:) = double(data.x.data(4,:));

N = length(x(1,:)); %length of data points there are


if cleaningON == true

    %First filter does a broadband filter
    f = 40; bw = 1.95;
    L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter **From WIKI: FIR filters express each output sample as a weighted sum of the last N input samples, where N is the order of the filter.
    B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter

    for i=1:4,
        r(i,:) = filtfilt(B, 1, x(i,:)); %r is filtered version of x (runs forward and backward using filtfilt -- if it goes in one direction, there will be a phase delay
    end

    %Next perform robust linear regression to find the best fit of tdTomato
    %onto mNeon and then substract
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
    warning('off','stats:statrobustfit:IterationLimit');
    numPeriods = 10;
    regressWin = round(window/(f*(1-0.5*bw))*numPeriods);
    tempz = z;
    z = tempocleanNorm(z,regressWin);
    warning('on','stats:statrobustfit:IterationLimit');

else %IF no signal cleaning - just apply a single filter around frequency of interest
    f = frequency;
    bw = bandwidth;
    L = ceil(3*(1/f*(1-0.5*bw))/dt+1); %length of time it will filter
    B = fir1(L, f*[1-0.5*bw 1+0.5*bw]*dt*2);%generates filter
    for i=1:4,
        z(i,:) = filtfilt(B, 1, x(i,:)); %z is filtered version of x (runs forward and backward using filtfilt -- if it goes in one direction, there will be a phase delay
    end
end

%Finally run through L and R mNeons and compare via cross-correlation.
%In addition compare to top 95% shuffled data and bottom 5% shuffled
%data

xCorrWin = floor(window/4); %Length of time (in samples) to perform cross correlations on
xCorrOvlp = 0;
xCorrStep = floor(window/4);

totalSS = 360; %Total amount of shift (goes half this left shifted (-) and half right shifted (+))
numSS = 16; %Total steps to get to shift. Will split half + and - around 0 degrees
xCorrShft = ((window/f)/(360/totalSS))/numSS;

%Note: Channel 3 (Right mNeon) gets shifted left (-) or right (+)
%Going through vector from first timepoint to end. take data in 1sec chunks. n is chunk number
n = 0;
for in=xCorrWin:xCorrStep:round(N-xCorrWin-xCorrOvlp-(xCorrShft*numSS))
        
    n=n+1;
    cor1(1+numSS,n) = xcorr(z(sigCHs(1),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(sigCHs(2),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),0,'normalized'); %Perform cross correlation between both adjusted and filtered mNeons
    corTD(1+numSS,n) = xcorr(z(ctrlCHs(1),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(ctrlCHs(2),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),0,'normalized'); %Perform cross correlation between left mNeon and right tdTomato - for investigating level of tdTomato contamination still present
   
    
    for j=1:numSS
        cor1(numSS-j+1,n) = xcorr(z(sigCHs(1),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(sigCHs(2),round(in-xCorrOvlp-(j*xCorrShft)):round(in+xCorrWin+xCorrOvlp-(j*xCorrShft))),0,'normalized'); %Perform cross correlation with a phase shift
        corTD(numSS-j+1,n) = xcorr(z(ctrlCHs(1),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(ctrlCHs(2),round(in-xCorrOvlp-(j*xCorrShft)):round(in+xCorrWin+xCorrOvlp-(j*xCorrShft))),0,'normalized'); %Perform cross correlation with a phase shift
    end
    
    for j=1:numSS
        cor1(j+1+numSS,n) = xcorr(z(sigCHs(1),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(sigCHs(2),round(in-xCorrOvlp+(j*xCorrShft)):round(in+xCorrWin+xCorrOvlp+(j*xCorrShft))),0,'normalized'); %Perform cross correlation with a phase shift
        corTD(j+1+numSS,n) = xcorr(z(ctrlCHs(1),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),z(ctrlCHs(2),round(in-xCorrOvlp+(j*xCorrShft)):round(in+xCorrWin+xCorrOvlp+(j*xCorrShft))),0,'normalized'); %Perform cross correlation with a phase shift
    end

end

cor1SHF = zeros(size(cor1,1), size(cor1,2), nshuffles);

if shuffleON
    fprintf("\nShuffle #: ");
    for shuffel_iter = 1:nshuffles
        fprintf("%d ",shuffel_iter);
        n=0;
        shfSignal = z(sigCHs(2),:); %Shuffled form of the 2nd signal channel
        shfSignal = shfSignal(randperm(length(shfSignal))); %Fully shuffle data
        for in=xCorrWin:xCorrStep:round(N-xCorrWin-xCorrOvlp-(xCorrShft*numSS))
            n=n+1;
            cor1SHF(1+numSS,n,shuffel_iter) = xcorr(z(sigCHs(1),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),shfSignal(1,in-xCorrOvlp:in+xCorrWin+xCorrOvlp),0,'normalized'); %Perform cross correlation between both adjusted and filtered mNeons

            for j=1:numSS
                cor1SHF(numSS-j+1,n,shuffel_iter) = xcorr(z(sigCHs(1),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),shfSignal(1,round(in-xCorrOvlp-(j*xCorrShft)):round(in+xCorrWin+xCorrOvlp-(j*xCorrShft))),0,'normalized'); %Perform cross correlation with a phase shift
            end

            for j=1:numSS
                cor1SHF(j+1+numSS,n,shuffel_iter) = xcorr(z(sigCHs(1),in-xCorrOvlp:in+xCorrWin+xCorrOvlp),shfSignal(1,round(in-xCorrOvlp+(j*xCorrShft)):round(in+xCorrWin+xCorrOvlp+(j*xCorrShft))),0,'normalized'); %Perform cross correlation with a phase shift
            end

        end
    end
    fprintf("...Done!\n");
end

cor1SHF = mean(cor1SHF,3);


clear numShftSteps window tmp regressWin R R10 nshuffles numPeriods N n in i L B bw f data dt




            