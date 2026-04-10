function trialSectionTimes = setTrialSections()

timestampNames = ["Trial Start","Dig Start","Outcome","Trial End","ITI End"];
trialSectionTimes(1,:) = {'Pre-Decision','Dig','Outcome','Post-Outcome','ITI Start','ITI End','Max-PreDec'}; % These are pre-defined trial sections

fprintf('\nDefining Trial Section Ranges - e.g. ''Pre-Decision'' period. Press any key to start.\n\n');

pause();

defaultValues = {[1,-1,2,0],[2,-2,2,2],[3,-2,3,2],[3,0,4,0],[4,0,4,15],[5,-15,5,0],3}; %These are values I've found to be optimal for my behavior - AJP

fprintf(['Default values are defined as: \n\n' ...
    'Pre-Decision: %s + %d to %s + %d\n' ...
    'Dig: %s + %d to %s + %d\n' ...
    'Outcome: %s + %d to %s + %d\n' ...
    'Post-Outcome: %s + %d to %s + %d\n' ...
    'ITI Start: %s + %d to %s + %d\n' ...
    'ITI End: %s + %d to %s + %d\n' ...
    'Maximum Pre-Decision: %d seonds\n'], ...
    timestampNames(defaultValues{1}(1)),defaultValues{1}(2),timestampNames(defaultValues{1}(3)),defaultValues{1}(4), ...
    timestampNames(defaultValues{2}(1)),defaultValues{2}(2),timestampNames(defaultValues{2}(3)),defaultValues{2}(4), ...
    timestampNames(defaultValues{3}(1)),defaultValues{3}(2),timestampNames(defaultValues{3}(3)),defaultValues{3}(4), ...
    timestampNames(defaultValues{4}(1)),defaultValues{4}(2),timestampNames(defaultValues{4}(3)),defaultValues{4}(4), ...
    timestampNames(defaultValues{5}(1)),defaultValues{5}(2),timestampNames(defaultValues{5}(3)),defaultValues{5}(4), ...
    timestampNames(defaultValues{6}(1)),defaultValues{6}(2),timestampNames(defaultValues{6}(3)),defaultValues{6}(4), ...
    defaultValues{7});

pendingInput = true;
while pendingInput
    userInput = input('Would you like to use the default values? (Y/N) ',"s");
    if strcmp(userInput,'N') || strcmp(userInput,'n')
        fprintf('Default values rejected - please enter new trial section ranges...\n')
        defaultValuesRejected = true;
        pendingInput = false;
    elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
        fprintf('Default values accepted!\n')
        defaultValuesRejected = false;
        pendingInput = false;
    else
        fprintf('Invalid input, please reenter values...\n');
    end
end

if defaultValuesRejected

    for sec = 1:6 %Pre-Dec through ITI End
        fprintf('\nDefining %s\n',trialSectionTimes{1,sec})
        pendingInput = true;
        while pendingInput
            timestamp1 = input(sprintf('%s starts in relation to which timestamp? (TStart = 1; Dig = 2; Outcome = 3; TEnd = 4; ITIEnd = 5) ',trialSectionTimes{1,sec}));
            offset1 = input(sprintf('%s starts how many seconds after %s? (use - to start before) ',trialSectionTimes{1,sec},timestampNames(timestamp1)));
            timestamp2 = input(sprintf('%s ends in relation to which timestamp? (TStart = 1; Dig = 2; Outcome = 3; TEnd = 4; ITIEnd = 5) ',trialSectionTimes{1,sec}));
            offset2 = input(sprintf('%s ends how many seconds after %s? (use - to end before) ',trialSectionTimes{1,sec},timestampNames(timestamp2)));

            userInput = input(sprintf('\n%s will be defined as %s + %d to %s + %d is this correct? (Y/N) ',trialSectionTimes{1,sec},timestampNames(timestamp1),offset1,timestampNames(timestamp2),offset2),"s");
            if strcmp(userInput,'N') || strcmp(userInput,'n')
                fprintf('Please reenter values...\n');
            elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
                pendingInput = false;
            else
                fprintf('Invalid input, please reenter values...\n');
            end
        end

        trialSectionTimes(2,sec) = {[timestamp1, offset1, timestamp2, offset2]};

    end

    trialSectionTimes(2,7) = {input('\nPlease select a maximum period for Pre-Decision (in seconds): ')};

else
    trialSectionTimes(2,:) = defaultValues;


end

end

