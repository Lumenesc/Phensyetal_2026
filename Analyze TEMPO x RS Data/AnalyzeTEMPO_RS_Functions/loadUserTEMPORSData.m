function cohortData = loadUserTEMPORSData()

%{
loadUserTEMPORSData v1 - AJP 11-2024

This is a console-user interface script that allows the user to load a
sigMatrix and timepoints matrix for each cohort they want to analyze


output: cohortData is a matrix consisting of 5 rows:

row 1: cohort 'name' provided by user
row 2: sigMatrix - loaded by user
row 3: timepoints - loaded by user and converted from seconds to 'samples'
row 4: sigMatrixCmplx - the converted sigMatrix into a complex vector
(explained in generateSigMatrixComplex function / documentation).
row 5: parameter metadata associated with the sigMatrixCmplx conversion

%}


%% Load user data - sigMatrix and timepoints / cohort

fprintf('Ok, Let us load your data first - you should have already generated a sigMatrix structure and timepoints dataset and saved these to your computer.\n');

pendingAnswer = true;
while pendingAnswer
    userInput = input('\nReady to load data? (Y/N): ',"s");
    if strcmp(userInput,'N') || strcmp(userInput,'n')
        fprintf('Please re-run script when ready to load data. Exiting now...\n');
        error('User cancelled data load.')
    elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
        pendingAnswer = false;
    else
        fprintf('Invalid Input\n');
    end

end

clear cohortData
loadCohorts = true; cohortCount = 0;
while loadCohorts

    cohortCount = cohortCount + 1;

    fprintf('\nLoading Cohort %d\n',cohortCount);
    userInput = input('Please enter a name for this cohort: ',"s");

    cohortData{cohortCount,1} = userInput;

    fprintf('Loading sigMatrix for %s...',cohortData{cohortCount,1});
    [file,path] = uigetfile({'*.mat'}, ...
        sprintf('Load sigMatrix Data File for %s\n',cohortData{cohortCount,1}));
    cohortData{cohortCount,2} = importdata([path,file]);

    fprintf('Loading timepoints for %s...',cohortData{cohortCount,1});
    [file,path] = uigetfile({'*.mat'}, ...
        sprintf('Load timepoints Data File for %s\n',cohortData{cohortCount,1}));
    cohortData{cohortCount,3} = importdata([path,file]);

    userInput = input('What was the step size used to build the sigMatrix (in ms)? (e.g. ''xCorrStep'')?');
    cohortData{cohortCount,3} = padITIEnd(cohortData{cohortCount,3},10); %Sometimes ITI End is missing on the last trial, this sets it to 10 seconds after the end
    cohortData{cohortCount,3} = convTimepointFS(cohortData{cohortCount,3},1000/userInput); %Converts the timepoints file from seconds to sigMatrix samples to appropriately index sigMatrix

    fprintf('sigMatrix and timepoints loaded for %s\n',cohortData{cohortCount,1});

    pendingAnswer = true;
    while pendingAnswer
        userInput = input('Would you like to load in another cohort? (Y,N)',"s");
        if strcmp(userInput,'N') || strcmp(userInput,'n')
            fprintf('Finished Loading Cohorts.\n');
            loadCohorts = false;
            pendingAnswer = false;
        elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end


end

%% Convert sigMatrix to Complex Vectors: sigMatrixCmplx

fprintf(['\nWe will now convert sigMatrix arrays from %d phase-shifted correlation values \n ' ...
    'to a complex vector consisting of the peak correlation value and phase for each sample\n'],size(cohortData{1,2}{1,3},1));
fprintf('\nFirst, we need to define the range of phase shifts performed.\n\n');

pendingAnswer = true;
while pendingAnswer
    MF_totalSS = input('What was the maximum degrees shifted? e.g.: ''totalSS''? ');
    MF_xCorrShft = input('How many degrees was each offset shifted? e.g. ''xCorrShft''? ');

    userInput = input(sprintf('\n sigMatrix calculated from %d to %d with %.02f increments. Is this correct? (Y/N)? ',-MF_totalSS,MF_totalSS,MF_xCorrShft),'s');
    if strcmp(userInput,'N') || strcmp(userInput,'n')
        fprintf('Please reenter parameters.\n');
    elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
        pendingAnswer = false;
    else
        fprintf('Invalid Input\n');
    end
end


MF_poIDX = [-MF_totalSS:MF_xCorrShft:MF_totalSS];

fprintf('\nNext let''s get the parameters for converting sigMatrix (Please check documentation for explanation of parameters):\n');
MF_hSW = input('Number of Samples for Horizontal Smoothing? (0 for no smoothing)? ');
MF_vSW = input('Number of Samples for Vertical Smoothing? (0 for no smoothing)? ');
thresh = input('What %tile threshold would you like to require for the peak correlation? (0 for no threshold, max 100)? ');

fprintf('\nComplex vectors will now be calculated for all sigMatrices provided - this may take a while.\n')

for cohortCount = 1:size(cohortData,1)
    [cohortData{cohortCount,4},cohortData{cohortCount,5}] = generateSigMatrixComplex(cohortData{cohortCount,1},cohortData{cohortCount,2},MF_poIDX,MF_hSW,MF_vSW,thresh);
end



end