function output = plotoscsynch(oscdata,oscChs,StateConstraints,IDX_SELECTION,IDX_RANGES,OMIT_CE,plotThreshold,removeZeros,plotPercent,performStats)

%%INPUTS: 




%Unpack oscdata
CText = oscdata{1};
stateArrayFull = oscdata{2};
%plotPhases =  oscdata{3};
plotSynch = oscdata{4};
IDX_ALL_INFO = oscdata{5};
IDX_ALL_CAT = oscdata{6};

figure

%Generate T - Table containing all combinations of the IDX inputs
T=combinations(IDX_RANGES{:});
if ~isempty(OMIT_CE)
    for i=1:size(T,1)
        for j=1:size(OMIT_CE,1)
            idx = find(IDX_SELECTION==OMIT_CE(j,1));
            if T{i,idx}==OMIT_CE(j,2)
                idx = find(IDX_SELECTION==3);
                T{i,idx} = 99;
            end
        end
    end
end
T=unique(T); T=standardizeMissing(T,99); %Remove duplicate rows and then make 99s NaN - note: unique does not work with NaN values, thus this step is required


%%%%%%%Setting Labels/Titles%%%%%%%

[PLOTTITLE,IDX_TITLES] = getIDXTitles(T,IDX_ALL_INFO,IDX_SELECTION,IDX_RANGES); %Plot title consists of all IDX selections that are of single length, IDX_TITLES contain those with more than one selection and will be displayed in the Legend
[SCTITLE,StateConstraintLabels] = getSCTitles(CText,StateConstraints); %State Constraint labels are pulled out of CText
[CHTITLE,THRESHTITLE,CHANNELNAMES] = getCHTitles(oscChs,plotThreshold);
[legendTitles] = getLegendTitles(IDX_TITLES,StateConstraintLabels,CHANNELNAMES,oscChs,T);

%Build a logic array to combine all state constraints
combStateIDX = ones(length(StateConstraints),size(stateArrayFull,2));
for i=1:length(StateConstraints)
    for j=1:length(StateConstraints{i})
        if StateConstraints{i}(j) ~= 0
            combStateIDX(i,:) = combStateIDX(i,:) & stateArrayFull(StateConstraints{i}(j),:);
        end
    end
end

output = cell.empty(0,3);
n=0;
for ch=1:length(oscChs)
    for con=1:length(StateConstraints)
        for i=1:size(T,1)
            n=n+1;

            %Build a logic array to combine all IDX constraints for current row of T
            values = rmmissing(T{i,:});
            combIDX = ones(1,size(IDX_ALL_CAT,2));
            for j=1:size(values,2)
                combIDX = combIDX & IDX_ALL_CAT(IDX_SELECTION(j),:)==IDX_ALL_INFO{2,3}{IDX_SELECTION(j)}(values(j));
            end

            %Apply the combIDX and combStateIDX to the plotSynch for the given channels
            data = getData(oscChs{ch},plotSynch,combIDX,combStateIDX(con,:));

            if removeZeros == true
                data = data(data~=0);
            elseif plotPercent==true
                data(data>0) = 1;
            end

            %Export data into output
            output{n,1} = legendTitles{n};
            output{n,2} = data';
            output{n,3} = combIDX & combStateIDX(con,:);

            %Plot the # of (or fraction if plotPrecent==true) timepoints exist for each theta bin as a polarplot            
            bar(n,mean(data))
            hold on
            errorbar(n,mean(data),std(data)/sqrt(length(data)),std(data)/sqrt(length(data)),'HandleVisibility','off')
            title(sprintf('%s when %s\n%s\n%s',CHTITLE,SCTITLE,PLOTTITLE,THRESHTITLE));
        end
    end
end

legend(legendTitles,'Location','southoutside')

if size(output,1) == 2 && performStats%if there are only two comparisons, perform a standard t test
    [~,p] = ttest2(output{1,2},output{2,2});
    fprintf('\nPaired T-Test between %s and %s: p = %.4f\n',output{1,1},output{2,1},p);
    annotation('textbox', [0.82, 0.15, 0.1, 0.1], 'String', sprintf("p=%.3f",p),'EdgeColor','none');
elseif performStats %Otherwise perform a one way ANOVA
    ANOVAdata = combCells(output(:,2)','double','post');
    [~,~,ANOVAmltpcmprsn]=runANOVA(ANOVAdata,'on');
    output{1,4} = ANOVAdata;
    output{1,5} = ANOVAmltpcmprsn;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUPPORT FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PLOTTITLE,IDX_TITLES] = getIDXTitles(T,IDX_ALL_INFO,IDX_SELECTION,IDX_RANGES)

for i=1:size(T,1)
    values = rmmissing(T{i,:});
    IDX_TITLES(i) = "";
    for j=1:size(values,2)
        if length(unique(T{:,j}))>1
            IDX_TITLES(i) = IDX_TITLES(i) + IDX_ALL_INFO{2,2}{IDX_SELECTION(j)}{values(j)} + " | ";
        end
    end
    if strlength(IDX_TITLES(i))>0
        IDX_TITLES(i) = extractBefore(IDX_TITLES(i),strlength(IDX_TITLES(i))-2);
    end
end

PLOTTITLE = "";
for i = 1:length(IDX_RANGES)
    if isscalar(IDX_RANGES{i})
        PLOTTITLE = PLOTTITLE + IDX_ALL_INFO{2,2}{IDX_SELECTION(i)}{IDX_RANGES{i}} + " | ";
    end
end
if strlength(PLOTTITLE)>0
    PLOTTITLE = "Trial Classification: " + extractBefore(PLOTTITLE,strlength(PLOTTITLE)-2); %Remove the last " | "
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SCTITLE,StateConstraintLabels] = getSCTitles(CText,StateConstraints)

%Get the titles for the State Constraints
StateConstraintLabels = cell.empty(0,length(StateConstraints));
for con=1:length(StateConstraints)
    StateConstraintLabels{con} = "";
    for cstate=StateConstraints{con}
        if cstate==0
            StateConstraintLabels{con}="Not Constrained";
        else
            for i=1:size(CText,2)-2
                if ~isempty(CText{cstate,i})
                    if StateConstraintLabels{con} == ""
                        StateConstraintLabels{con} = CText{cstate,i};
                    else
                        StateConstraintLabels{con} = StateConstraintLabels{con} + " & " + CText{cstate,i};
                    end
                end
            end
        end
        if  StateConstraintLabels{con}==""
            StateConstraintLabels{con}="| No Synch";
        end
    end
end
if isempty(StateConstraints)
    SCTITLE = 'No Constraints';
elseif isscalar(StateConstraintLabels)
    SCTITLE = StateConstraintLabels{1};
else
    SCTITLE = "Constrained";
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CHTITLE,THRESHTITLE,CHANNELNAMES] = getCHTitles(oscChs,plotThreshold)

CHANNELNAMES = string.empty(0,length(oscChs));
for i=1:length(oscChs)
    if strcmp(oscChs{i},'wH')
        CHANNELNAMES{i} = 'wH-PVxMD';
    elseif strcmp(oscChs{i},'xH')
        CHANNELNAMES{i} = 'xH-PVxMD';
    elseif strcmp(oscChs{i},'PV')
        CHANNELNAMES{i} = 'PVxPV';
    elseif strcmp(oscChs{i},'MD')
        CHANNELNAMES{i} = 'MDxMD';
    elseif strcmp(oscChs{i},'wH-PD')
        CHANNELNAMES{i} = 'wH-PVxMD Phase Difference';
    elseif strcmp(oscChs{i},'xH-PD')
        CHANNELNAMES{i} = 'xH-PVxMD Phase Difference';
    end
end
if isscalar(CHANNELNAMES)
    CHTITLE = CHANNELNAMES{1};
else
    CHTITLE = 'Synchrony';
end

if plotThreshold==true
    THRESHTITLE = "Threshold: N/A";
else
    THRESHTITLE = "Threshold: N/A";
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function legendTitles = getLegendTitles(IDX_TITLES,SCLABELS,CHANNELNAMES,oscChs,T)

legendTitles = cell.empty(0,size(T,1)); n=0;
for ch=1:length(oscChs)
    for con=1:length(SCLABELS)
        for i=1:size(T,1)
            n=n+1;
            if ~isscalar(CHANNELNAMES) && ~isscalar(SCLABELS) && strlength(IDX_TITLES(i))>0
                legendTitles{1,n} = CHANNELNAMES{ch} + " : when " + SCLABELS{con} + " : " + IDX_TITLES(i);
            elseif ~isscalar(CHANNELNAMES) && strlength(IDX_TITLES(i))>0
                legendTitles{1,n} = CHANNELNAMES{ch} + " : " + IDX_TITLES(i);
            elseif ~isscalar(SCLABELS) && strlength(IDX_TITLES(i))>0
                legendTitles{1,n} = "when " + SCLABELS{con} + " : " + IDX_TITLES(i);
            elseif ~isscalar(CHANNELNAMES) && ~isscalar(SCLABELS)
                legendTitles{1,n} = CHANNELNAMES{ch} + " : when " + SCLABELS{con};
            elseif ~isscalar(CHANNELNAMES)
                legendTitles{1,n} = CHANNELNAMES{ch};
            elseif  ~isscalar(SCLABELS)
                legendTitles{1,n} = "when " + SCLABELS{con};
            elseif strlength(IDX_TITLES(i))>0
                legendTitles{1,n} = IDX_TITLES(i);
            else
                legendTitles{1,n} = "";
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = getData(currentCH,plotSynch,combIDX,combStateIDX)

%CHANNELS - MANUALLY DEFINED, POTENTIALL ALLOW FOR MORE FLEXIBILITY IN THE FUTURE

%From [[PVL,MDL,~,~,~];[PVL,PVR,~,~,~];[PVL,MDR,~,~,~];[PVR,MDL,~,~,~];[MDL,MDR,~,~,~];[PVR,MDR,~,~,~]]
wHLeft = 1;
PV = 2;
xHLeftPV = 3;
xHRightPV = 4;
MD = 5;
wHRight = 6;

%Select the appropriate rows from plotSynch
if strcmp(currentCH,'xH')
        data = [plotSynch(xHLeftPV,combIDX & combStateIDX) plotSynch(xHRightPV,combIDX & combStateIDX)];
elseif strcmp(currentCH,'wH')
        data = [plotSynch(wHLeft,(combIDX & combStateIDX)) plotSynch(wHRight,(combIDX & combStateIDX))];
elseif strcmp(currentCH,'PV')
        data = plotSynch(PV,combIDX & combStateIDX);
elseif strcmp(currentCH,'MD')
        data = plotSynch(MD,combIDX & combStateIDX);
elseif strcmp(currentCH,'xH-PD')
        data = wrapTo2Pi(abs(plotSynch(xHLeftPV,combIDX & combStateIDX)- plotSynch(xHRightPV,combIDX & combStateIDX)));
elseif strcmp(currentCH,'wH-PD')
        data = wrapTo2Pi(abs(plotSynch(wHLeft,combIDX & combStateIDX)- plotSynch(wHRight,combIDX & combStateIDX)));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

