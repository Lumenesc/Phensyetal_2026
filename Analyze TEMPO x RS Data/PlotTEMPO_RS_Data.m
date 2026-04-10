
%{

PlotTEMPO_RS_Data v1 

Requires: AnalyzeTEMPO_RS_Data to be run prior.

This script serves as the master script to handle plotting of TEMPO x RS
Synchrony data. The user will be able to first select which trial
grouping/classification they would prefer to examine (or a combination) and
then plot the synchrony data within the selected group. By default each
cohort will be plotted separately. 

Note: the trial groupings are defined as indices and thus are referred to 
as IDXs.


The following types of plots are currently supported:

1. Bar Plots: Average Synchrony of Trial Sections within IDX(s)




Note: the script allows for the exporting of the data in the plots in case
the user would like to generate their own graphs.

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
%Setup / clearup
close all 
for i=1:15
    for j=1:i
        fprintf('*');
    end
    fprintf('\n');
end
fprintf('Executing PlotTEMPO_RS_Functions\n\n');
fprintf(['This script will help you graphically examine the TEMPO synchrony data calculated in AnalyzeTEMPO_RS_Functions.\n' ...
    'Since Rule Shift data can be organized into many different trial groupings, this script allows you to examine a number\n' ...
    'of different pre-defined groupings or combinations (e.g. correct/error + conflict/nonconflict trials). You will also be\n' ...
    'prompted to export the data in the graphs after plotting is complete.\n'])
fprintf('\nPress any key to get started.\n');
pause();
clear i j

addpath('distFig v5.5')
addpath('AnalyzeTEMPO_RS_Functions')

if ~exist("synchronyData",'var')
    error('Synchrony Data does not exist - please run AnalyzeTEMPO_RS_Data first.')
end
if ~exist("behavioralIDXData",'var')
    error('Trial IDX Grouping Data does not exist - please run AnalyzeTEMPO_RS_Data first.')
end



%% Choose which Analysis Group(s) to plot

%%%%%%%%NEED TO ADD IN FUNCTIONALITY FOR MULTIPLE PHASE RANGES%%%%%%%%%%%%

if exist('uselect_analysisGroup','var')
    pendingAnswer = true;
    while pendingAnswer
        fprintf('\nCurrent Groups being Analyzed:\n')
        for analysisGroup=uselect_analysisGroup
            fprintf('%d. %s\n',analysisGroup,synchronyData{1,analysisGroup})
        end
        fprintf('\n');
        userInput = input('Keep current selection? (Y/N): ',"s");
        if strcmp(userInput,'Y') || strcmp(userInput,'y')
            fprintf('Continuing with current analysis group(s)...\n\n');
            pendingAnswer = false;
        elseif strcmp(userInput,'N') || strcmp(userInput,'n')
            uselect_analysisGroup = selectAnalysisGroups(synchronyData);
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
else
    uselect_analysisGroup = selectAnalysisGroups(synchronyData);
end

if exist('uselect_phaseset','var')
    pendingAnswer = true;
    while pendingAnswer
        fprintf('\nThe following phase ranges are selected for plotting:\n')
        for phrng=uselect_phaseset
            fprintf('%d. %s\n',phrng,synchronyData{2,1}{1,phrng})
        end
        fprintf('\n');
        userInput = input('Keep current selection? (Y/N): ',"s");
        if strcmp(userInput,'Y') || strcmp(userInput,'y')
            fprintf('Continuing with current phase range(s)...\n\n');
            pendingAnswer = false;
        elseif strcmp(userInput,'N') || strcmp(userInput,'n')
            uselect_phaseset = selectPhaseRanges(synchronyData);
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
else
    uselect_phaseset = selectPhaseRanges(synchronyData);
end

clear phrng pendingAnswer userInput

%% Choose IDX Categories


if exist('uselect_idx','var')
    pendingAnswer = true;
    while pendingAnswer
        fprintf('\nCurrent category selection:\n')
        for IDX=uselect_idx
            fprintf('%d. %s\n',IDX,behavioralIDXData{2,1}{1,IDX})
        end
        fprintf('\n');
        userInput = input('Keep current selection? (Y/N): ',"s");
        if strcmp(userInput,'Y') || strcmp(userInput,'y')
            fprintf('Continuing with current category selection...\n\n');
            pendingAnswer = false;
        elseif strcmp(userInput,'N') || strcmp(userInput,'n')
            uselect_idx = getIDXSelection(behavioralIDXData);
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
else
    uselect_idx = getIDXSelection(behavioralIDXData);
end

clear IDX pendingAnswer userInput
%% Choose Task Phase and Trial Sections

if exist('uselect_taskphase','var') && exist('uselect_trsections','var')
    pendingAnswer = true;
    while pendingAnswer
        fprintf('\nCurrent task selections:\n\n')
        fprintf('Task Phase(s):\n');
            for taskphase=uselect_taskphase
                if taskphase==1
                    fprintf('1. Initial Association\n')
                elseif taskphase==2
                    fprintf('2. Rule Shift\n')
                end
            end
        fprintf('\nTrial Sections:\n')
        for trsec=uselect_trsections
            fprintf('%d. %s\n',trsec,behavioralIDXData{3,1}{trsec})
        end
        fprintf('\n');
        userInput = input('Keep current selection? (Y/N): ',"s");
        if strcmp(userInput,'Y') || strcmp(userInput,'y')
            fprintf('Continuing with current task selections...\n\n');
            pendingAnswer = false;
        elseif strcmp(userInput,'N') || strcmp(userInput,'n')
            [uselect_taskphase,uselect_trsections] = getTaskSelections(behavioralIDXData);
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
else
    [uselect_taskphase,uselect_trsections] = getTaskSelections(behavioralIDXData);
end

clear trsec pendingAnswer userInput

%% Plot Formatting

plotLevels = loadPlotLevels(synchronyData,behavioralIDXData,uselect_analysisGroup,uselect_phaseset,uselect_taskphase,uselect_trsections,uselect_idx);

if exist('uselect_plotFormat','var')
    pendingAnswer = true;
    while pendingAnswer
        fprintf('\nCurrent plot format selection:\n')
        fprintf('\n%s will be plotted for each %s. The rest of the categories wil be organized into individual plots.\n',plotLevels{1,uselect_plotFormat(1)},plotLevels{1,uselect_plotFormat(2)});
        userInput = input('Keep current selection? (Y/N): ',"s");
        if strcmp(userInput,'Y') || strcmp(userInput,'y')
            fprintf('Continuing with current plot formatting...\n\n');
            pendingAnswer = false;
        elseif strcmp(userInput,'N') || strcmp(userInput,'n')
            uselect_plotFormat = getPlotFormat(plotLevels);
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
else
    uselect_plotFormat = getPlotFormat(plotLevels);
end

clear trsec pendingAnswer userInput


%% Do the plotting - UNDER CONSTRUCTION!!!
close all

ybound = [-2 2]; yTitle = 'Synchrony';

plotIDX = [1:size(plotLevels,2)];
plotIDX = [plotIDX(plotIDX~=uselect_plotFormat(1) & plotIDX~=uselect_plotFormat(2)),uselect_plotFormat(1),uselect_plotFormat(2)]; %Remove the 'base' and 'group' selections

plotLevelsSrt = plotLevels(:,plotIDX);
T = combinations(plotLevels{3,plotIDX});
T.Properties.VariableNames = plotLevels(1,plotIDX);

clear output
for trow = 1:size(T,1)
    [IDX_set,IDX_names,IDX_ranges] = buildIDXSet([T{trow,'Analysis Group'}],behavioralIDXData(2,:),uselect_idx);
    IDX_set = combineIDXSets(IDX_set,IDX_ranges);

    for k = 1:size(synchronyData{2,T{trow,'Analysis Group'}}{2,T{trow,'Phase Range'}},1)
        IDX_Mask = IDX_set{T{trow,'Trial Grouping (IDX)'},T{trow,'Task Phase'}}{k};
        output(trow,k) = mean(synchronyData{2,T{trow,'Analysis Group'}}{2,T{trow,'Phase Range'}}{k,T{trow,'Task Phase'}}(IDX_Mask,T{trow,'Trial Section'}),'omitnan');
    end
end

close all
clear plotdata
T2 = unique(T(:,1:3));
for i=1:size(unique(T(:,1:3)),1)
    plotdata{1,i} = sprintf('%s - %s - %s', ...
        plotLevelsSrt{2,1}{T2{i,1}}, ...
        plotLevelsSrt{2,2}{T2{i,2}}, ...
        plotLevelsSrt{2,3}{T2{i,3}});

    n=0;
    for group = unique(T{:,4})'
        n=n+1;
        plotdata{2,i}{1,n} = plotLevelsSrt{2,4}{group};
        plotdata{2,i}{2,n} = output(T.(1)==T2{i,1} & T.(2)==T2{i,2} & T.(3)==T2{i,3} & T.(4)==group,:)';
    end

    plotgraph_bar(plotdata{2,i}(2,:),plotLevelsSrt{2,5}(plotLevelsSrt{3,5}),ybound,yTitle,plotdata{2,i}(1,:),[],plotdata{1,i},true)
end
distFig

clear i k n trow group output T T2 ybound yTitle IDX_Mask IDX_names IDX_ranges IDX_set plotIDX plotLevelsSrt

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZATION FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotLevels = loadPlotLevels(synchronyData,behavioralIDXData,uselect_analysisGroup,uselect_phaseset,uselect_taskphase,uselect_trsections,uselect_idx)

plotLevels = cell(0,5);
plotLevels(1,:) = {'Analysis Group','Phase Range','Task Phase','Trial Section','Trial Grouping (IDX)'};
plotLevels{2,1} = synchronyData(1,:);
plotLevels{2,2} = synchronyData{2,1}(1,:);
plotLevels{2,3} = {'IA','RS'};
plotLevels{2,4} = behavioralIDXData{3,1};

n=0;
for i=uselect_idx
    n=n+1;
    IDX_names{1,n} = behavioralIDXData{2,1}{3,i};
end
IDX_names = combineIDXNames(IDX_names); %Collate the separate IDX tags
plotLevels{2,5} = IDX_names;


plotLevels{3,1} = uselect_analysisGroup;
plotLevels{3,2} = uselect_phaseset;
plotLevels{3,3} = uselect_taskphase;
plotLevels{3,4} = uselect_trsections;
plotLevels{3,5} = 1:length(IDX_names);


end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  USER DATA COLLECTION FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uselect_analysisGroup = selectAnalysisGroups(synchronyData)
fprintf('\nGroups associated with analyzed TEMPO Synchrony dataset:\n')
for analysisGroup = 1:size(synchronyData,2)
    fprintf('%d. %s\n',analysisGroup,synchronyData{1,analysisGroup});
end

pendingAnswer = true;
while pendingAnswer
    userInput = (split(input('\nWhich group(s) would you like to plot? (Separate groups with a comma or space). ','s'),{',',' '}));
    if isempty(userInput)
        fprintf('No selection made, please try again.\n\n')
    elseif length(userInput) > size(synchronyData{2,1},2)
        fprintf('\nToo many entries, please try again.\n')
    else
        uselect_analysisGroup = double.empty(0,length(userInput));
        for i=1:length(userInput)
            uselect_analysisGroup(i) = str2double(cell2mat(userInput(i)));
        end
        if any(isnan(uselect_analysisGroup))
            fprintf('\nAt least one entry is invalid, please try again.\n')
        elseif length(uselect_analysisGroup)~=length(unique(uselect_analysisGroup))
            fprintf('\nCannot accept redundant entries, please try again.\n')
        else
            fprintf('Phase Range(s) Selected:\n');
            for analysisGroup=uselect_analysisGroup
                fprintf('%d. %s\n',analysisGroup,synchronyData{1,analysisGroup})
            end
            userInput = input('Is this correct? (Y/N) ','s');
            if strcmp(userInput,'N') || strcmp(userInput,'n')
                fprintf('Please reenter choice.\n');
            elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
                pendingAnswer = false;
            else
                fprintf('Invalid Input\n');
            end
        end
    end

    
    if isempty(userInput)
        fprintf('No selection made, please try again.\n\n')
    elseif any(~ismember(userInput,synchronyData(1,:)))
        fprintf('At least one selection does not match the groups analyzed, please try again.\n\n')
    elseif length(userInput)~=length(unique(userInput))
        fprintf('\nCannot accept redundant entries, please try again.\n')
    else
        [~,uselect_analysisGroup] = ismember(userInput,synchronyData(1,:));
        uselect_analysisGroup = uselect_analysisGroup';
        fprintf('Group(s) to plot:\n');
        for analysisGroup=uselect_analysisGroup
            fprintf('%d. %s\n',analysisGroup,synchronyData{1,analysisGroup})
        end
        userInput = input('Is this correct? (Y/N) ','s');
        if strcmp(userInput,'N') || strcmp(userInput,'n')
            fprintf('Please reenter choice.\n');
        elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
end
end

function uselect_phaseset = selectPhaseRanges(synchronyData)
fprintf('\nPhase ranges associated with analyzed TEMPO Synchrony dataset:\n')
for phrng = 1:size(synchronyData{2,1},2)
    fprintf('%d. %s\n',phrng,synchronyData{2,1}{1,phrng});
end

pendingAnswer = true;
while pendingAnswer
    userInput = (split(input('\nWhich group(s) would you like to plot? (Separate groups with a comma or space). ','s'),{',',' '}));
    if isempty(userInput)
        fprintf('No selection made, please try again.\n\n')
    elseif length(userInput) > size(synchronyData{2,1},2)
        fprintf('\nToo many entries, please try again.\n')
    else
        uselect_phaseset = double.empty(0,length(userInput));
        for i=1:length(userInput)
            uselect_phaseset(i) = str2double(cell2mat(userInput(i)));
        end
        if any(isnan(uselect_phaseset))
            fprintf('\nAt least one entry is invalid, please try again.\n')
        elseif length(uselect_phaseset)~=length(unique(uselect_phaseset))
            fprintf('\nCannot accept redundant entries, please try again.\n')
        else
            fprintf('Phase Range(s) Selected:\n');
            for phrng=uselect_phaseset
                fprintf('%d. %s\n',phrng,synchronyData{2,1}{1,phrng})
            end
            userInput = input('Is this correct? (Y/N) ','s');
            if strcmp(userInput,'N') || strcmp(userInput,'n')
                fprintf('Please reenter choice.\n');
            elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
                pendingAnswer = false;
            else
                fprintf('Invalid Input\n');
            end
        end
    end
end
end

function uselect_idx = getIDXSelection(behavioralIDXData)
fprintf('\nThe current trial groupings are as follows (see documentation for description):\n\n');
for IDX=1:size(behavioralIDXData{2,1},2)
    fprintf('%02d. %s\n',IDX,behavioralIDXData{2,1}{1,IDX})
end

pendingAnswer = true;
while pendingAnswer
    userInput = (split(input('\nWhich categorie(s) would you like to plot (up to 3 combinations)? (Separate categories with a comma or space). ','s'),{',',' '}));
    if isempty(userInput)
        fprintf('No selection made, please try again.\n\n')
    elseif length(userInput) > 3
        fprintf('\nToo many combinations - maximum is 3, please try again.\n')
    else
        uselect_idx = double.empty(0,length(userInput));
        for i=1:length(userInput)
            uselect_idx(i) = str2double(cell2mat(userInput(i)));
        end
        if any(isnan(uselect_idx))
            fprintf('\nAt least one entry is invalid, please try again.\n')
        elseif any(uselect_idx > size(behavioralIDXData{2,1},2))
            fprintf('\nAt least one entry is invalid, please try again.\n')
        elseif length(uselect_idx)~=length(unique(uselect_idx))
            fprintf('\nCannot accept redundant entries, please try again.\n')
        else
            fprintf('Category Selection(s):\n');
            for IDX=uselect_idx
                fprintf('%d. %s\n',IDX,behavioralIDXData{2,1}{1,IDX})
            end
            userInput = input('Is this correct? (Y/N) ','s');
            if strcmp(userInput,'N') || strcmp(userInput,'n')
                fprintf('Please reenter choice.\n');
            elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
                pendingAnswer = false;
            else
                fprintf('Invalid Input\n');
            end
        end
    end
end
end

function [uselect_taskphase,uselect_trsections] = getTaskSelections(behavioralIDXData)

pendingAnswer = true;
while pendingAnswer

    fprintf('The task is broken into two phases:\n\n')
    fprintf(['01. Initial Association\n' ...
        '02. Rule Shift']);

     userInput = (split(input('\nWhat task phase(s) do you want to plot? (Separate categories with a comma or space). ','s'),{',',' '}));
    if isempty(userInput)
        fprintf('No selection made, please try again.\n\n')
    else
        uselect_taskphase = double.empty(0,length(userInput));
        for i=1:length(userInput)
            uselect_taskphase(i) = str2double(cell2mat(userInput(i)));
        end
        if any(isnan(uselect_taskphase))
            fprintf('\nAt least one entry is invalid, please try again.\n')
        elseif any(uselect_taskphase > size(behavioralIDXData{3,1},2))
            fprintf('\nAt least one entry is invalid, please try again.\n')
        elseif length(uselect_taskphase)~=length(unique(uselect_taskphase))
            fprintf('\nCannot accept redundant entries, please try again.\n')
        else
            fprintf('Task Phase(s) Selected:\n');
            for taskphase=uselect_taskphase
                if taskphase==1
                    fprintf('1. Initial Association\n')
                elseif taskphase==2
                    fprintf('2. Rule Shift\n')
                end
            end
            userInput = input('Is this correct? (Y/N) ','s');
            if strcmp(userInput,'N') || strcmp(userInput,'n')
                fprintf('Please reenter choice.\n');
            elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
                pendingAnswer = false;
            else
                fprintf('Invalid Input\n');
            end
        end
    end

    if strcmp(userInput,'IA') || strcmp(userInput,'ia')
        uselect_taskphase = 1;
        pendingAnswer=false;
    elseif strcmp(userInput,'RS') || strcmp(userInput,'rs')
        uselect_taskphase = 2;
        pendingAnswer=false;
    else
        fprintf('Invalid Input\n');
    end
end

pendingAnswer = true;
while pendingAnswer
    fprintf('\nTrial sections available for plotting are as follows:\n')
    for tasksection = 1:length(behavioralIDXData{3,1})
        fprintf('%02d. %s\n',tasksection,behavioralIDXData{3,1}{tasksection})
    end
    userInput = (split(input('\nWhich trial section(s) would you like to plot? (Separate categories with a comma or space). ','s'),{',',' '}));
    if isempty(userInput)
        fprintf('No selection made, please try again.\n\n')
    else
        uselect_trsections = double.empty(0,length(userInput));
        for i=1:length(userInput)
            uselect_trsections(i) = str2double(cell2mat(userInput(i)));
        end
        if any(isnan(uselect_trsections))
            fprintf('\nAt least one entry is invalid, please try again.\n')
        elseif any(uselect_trsections > size(behavioralIDXData{3,1},2))
            fprintf('\nAt least one entry is invalid, please try again.\n')
        elseif length(uselect_trsections)~=length(unique(uselect_trsections))
            fprintf('\nCannot accept redundant entries, please try again.\n')
        else
            fprintf('Trial Sections Selected:\n');
            for trsec=uselect_trsections
                fprintf('%d. %s\n',trsec,behavioralIDXData{3,1}{trsec})
            end
            userInput = input('Is this correct? (Y/N) ','s');
            if strcmp(userInput,'N') || strcmp(userInput,'n')
                fprintf('Please reenter choice.\n');
            elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
                pendingAnswer = false;
            else
                fprintf('Invalid Input\n');
            end
        end
    end
end

end

function uselect_plotFormat = getPlotFormat(plotLevels)

fprintf(['\nWe will now decide how to plot the data. There are %d different levels to plot, however only 2\n' ...
    'levels can be organized into a single plot - a ''base'' and a ''group'' level plot. The rest of the\n' ...
    'levels will be organized into separate plots.\n'],size(plotLevels,2));
fprintf('\nPress any key when ready to set plot parameters.\n');
pause();

fprintf('\nThe ''base'' plot is the lowest plot level and will be plotted for each ''group''. The options for either selection are:\n\n')
for option=1:size(plotLevels,2)
    fprintf('%02d. %s\n',option,plotLevels{1,option})
end

pendingAnswer = true;
while pendingAnswer
    uselect_plotFormat(1) = input('\nWhich category would you like to use for the ''base'' plot? ');
    uselect_plotFormat(2) = input('Which category would you like to use for the ''group'' plot? ');

    if uselect_plotFormat(1) == uselect_plotFormat(2)
        fprintf('The ''base'' and ''group'' selections cannot be the same, please try again.\n');
    elseif uselect_plotFormat(1) > size(plotLevels,2) || uselect_plotFormat(1) < 1 || uselect_plotFormat(2) > size(plotLevels,2) || uselect_plotFormat(2) < 1
        fprintf('\nAt least one entry is invalid, please try again.\n')
    else
        fprintf('\n%s will be plotted for each %s. The rest of the categories wil be organized into individual plots.\n',plotLevels{1,uselect_plotFormat(1)},plotLevels{1,uselect_plotFormat(2)});
        userInput = input('Is this correct? (Y/N) ','s');
        if strcmp(userInput,'N') || strcmp(userInput,'n')
            fprintf('Please reenter choice.\n');
        elseif strcmp(userInput,'Y') || strcmp(userInput,'y')
            pendingAnswer = false;
        else
            fprintf('Invalid Input\n');
        end
    end
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  SUPPORT FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IDX_set,IDX_names,IDX_ranges] = buildIDXSet(cohortSet,IDX_set_cohort,idx_subset)

IA = 1; RS = 2;
IDX_set = {}; IDX_names={}; IDX_ranges={};
for i=1:size(idx_subset,2)

    IDX_sub_IA = cell(0,1);
    IDX_sub_RS = cell(0,1);

    for cohort=cohortSet
        IDX_sub_IA = [IDX_sub_IA;IDX_set_cohort{1,cohort}{2,idx_subset(i)}(:,IA)];
        IDX_sub_RS = [IDX_sub_RS;IDX_set_cohort{1,cohort}{2,idx_subset(i)}(:,RS)];
    end
    IDX_set{i,IA} = IDX_sub_IA;
    IDX_set{i,RS} = IDX_sub_RS;

    IDX_names{1,i} = IDX_set_cohort{1,cohort}{3,idx_subset(i)};
    IDX_ranges{1,i} = IDX_set_cohort{1,cohort}{4,idx_subset(i)};

end

IDX_names = combineIDXNames(IDX_names); %Collate the separate IDX tags

end

function IDX_names_combined = combineIDXNames(IDX_name_set)

n=0;

if length(IDX_name_set) > 3
    error('combine IDX Names cannot handle a name set of that length');
elseif isempty(IDX_name_set)
    IDX_names_combined = {'No sub-index'};
elseif length(IDX_name_set) == 1
    IDX_names_combined = IDX_name_set{1,1};
elseif length(IDX_name_set) == 2
    for i=1:length(IDX_name_set{1})
        for j=1:length(IDX_name_set{2})
            n=n+1;
            IDX_names_combined{n} = sprintf('%s + %s',IDX_name_set{1}{i},IDX_name_set{2}{j});
        end
    end
elseif length(IDX_name_set) == 3
    for i=1:length(IDX_name_set{1})
        for j=1:length(IDX_name_set{2})
            for m=1:length(IDX_name_set{3})
                n=n+1;
                IDX_names_combined{n} = sprintf('%s + %s + %s',IDX_name_set{1}{i},IDX_name_set{2}{j},IDX_name_set{3}{m});
            end
        end
    end
end

end

function IDX_set_out = combineIDXSets(IDX_set,IDX_ranges)

if size(IDX_ranges,2) > 3
    error('Too many IDX selections - cannot handle more than 3');
end

IA=1; RS=2;

for uselect_taskphase = IA:RS
    n=0;
    for i=IDX_ranges{1,1}
        if size(IDX_ranges,2) > 1
            for j=IDX_ranges{1,2}
                if size(IDX_ranges,2) > 2
                    for h=IDX_ranges{1,3}
                        n=n+1;
                        for k=1:size(IDX_set{1,uselect_taskphase},1)
                            IDX_set_out{n,uselect_taskphase}{k,1} = IDX_set{1,uselect_taskphase}{k,1}==i & IDX_set{2,uselect_taskphase}{k,1}==j & IDX_set{3,uselect_taskphase}{k,1}==h;
                        end
                    end
                else
                    n=n+1;
                    for k=1:size(IDX_set{1,uselect_taskphase},1)
                        IDX_set_out{n,uselect_taskphase}{k,1} = IDX_set{1,uselect_taskphase}{k,1}==i & IDX_set{2,uselect_taskphase}{k,1}==j;
                    end
                end
            end
        else
            n=n+1;
            for k=1:size(IDX_set{1,uselect_taskphase},1)
                IDX_set_out{n,uselect_taskphase}{k,1} = IDX_set{1,uselect_taskphase}{k,1}==i;
            end
        end
    end

end

end
