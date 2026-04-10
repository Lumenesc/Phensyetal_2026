%% Builds the 'masterx' which filters, cleans, and concatenates all animal's data.x.data in a folder

%First - select folder to populate each animal's session length (in
%samples)
files = dir('C:\Users\AJPhe\Box\Postdoctoral Work\Gamma Oscillation Project\TEMPO Data\H2C\Day 1');
files(ismember( {files.name}, {'.', '..'})) = [];  %remove . and ..
dirFlags = [files.isdir]; % Get a logical vector that tells which is a directory.
subFolders = files(dirFlags); % Extract only those that are directories.

%Collect session length for each subdirectory folder ##NOTE: All subdirectories must be SEV container folders
masterx = double.empty(2,0);
masterxREF = double.empty(2,0);
for k = 1 : length(subFolders)
    [data] = SEV2mat(strcat(subFolders(k).folder,'\',subFolders(k).name));
    fprintf('Collecting data on sub folder #%d/%d = %s...\n', k, length(subFolders), subFolders(k).name);



    frequency = 40; % frequency to analyze (in Hz)
    bandwidth = (50-30) / frequency; % filter around frequency using this
    window = round(data.x.fs); %fs is sampling frequencies. window is the number of samples in each second (window is 1sec)
    dt = 1/data.x.fs; %timestep

    clear x z r
    %convert datatype to double with sev function
    x(1,:) = double(data.x.data(1,:)); %Left ACE-mNeon
    x(2,:) = double(data.x.data(2,:)); %Left Varnam
    x(3,:) = double(data.x.data(3,:)); %Right ACE-mNeon
    x(4,:) = double(data.x.data(4,:)); %Right Varnam
    x(5,:) = double(data.x.data(5,:)); %Left cyOFP
    x(6,:) = double(data.x.data(6,:)); %Right cyOFP

    N = length(x(1,:)); %length of data points there are


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

    masterx = [masterx z(1:4,:)];
    masterxREF = [masterxREF z(5:6,:)];
end
clear files dirFlags subFolders k