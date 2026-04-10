function output = plotoscphaseheatmap(oscdata,oscChs,oscChLabels,StateConstraints,IDX_SELECTION,IDX_RANGES,OMIT_CE,plotPercent,plotThreshold,phaseBinIncrement,phaseBinWidth,histogramON,performWatson)

%%INPUTS:

% if length(oscChs) ~= 2
%     error("plotoscphaseheatmap requires two oscChs to use as heatmap axises.");
% end

MAX_COUNT = 0;

%Unpack oscdata
CText = oscdata{1};
stateArrayFull = oscdata{2};
plotPhases = wrapTo2Pi(oscdata{3});
%plotSynch = oscdata{4};
IDX_ALL_INFO = oscdata{5};
IDX_ALL_CAT = oscdata{6};

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
[CHTITLE,THRESHTITLE,CHANNELNAMES] = getCHTitles(oscChs{end},plotThreshold);
%[legendTitles] = getLegendTitles(IDX_TITLES,StateConstraintLabels,CHANNELNAMES,oscChs,T);

%Build a logic array to combine all state constraints
combStateIDX = ones(length(StateConstraints),size(stateArrayFull,2));
for i=1:length(StateConstraints)
    for j=1:length(StateConstraints{i})
        if StateConstraints{i}(j) ~= 0
            combStateIDX(i,:) = combStateIDX(i,:) & stateArrayFull(StateConstraints{i}(j),:);
        end
    end
end

output = cell.empty(0,2);
n=0;
for con=1:length(StateConstraints)
    for i=1:size(T,1)
        n=n+1;
        %Build a logic array to combine all IDX constraints for current row of T
        values = rmmissing(T{i,:});
        combIDX = ones(1,size(IDX_ALL_CAT,2));
        for j=1:size(values,2)
            combIDX = combIDX & IDX_ALL_CAT(IDX_SELECTION(j),:)==IDX_ALL_INFO{2,3}{IDX_SELECTION(j)}(values(j));
        end
        for ch=1:length(oscChs)
            for subch=1:length(oscChs{ch})
                %Apply the combIDX and combStateIDX to the plotPhases for the given channels
                data{ch,subch} = getData(oscChs{ch}{subch},plotPhases,combIDX,combStateIDX(con,:));
            end
        end
        f(n) = figure;
        sgtitle(sprintf('%s when %s\n%s: %s\n%s',CHTITLE,SCTITLE,PLOTTITLE,IDX_TITLES{i},THRESHTITLE));
        for ch=1:length(oscChs)-1
            %Calculate the frequency of timepoints within each provided phase range (based on phaseBinIncrement and phaseBinRange)
            thetas = deg2rad([135:phaseBinIncrement:225]); %thetas = circshift(thetas,round(length(thetas)/2));
            thetas2 = deg2rad([-180:phaseBinIncrement:180]); %thetas2 = circshift(thetas2,round(length(thetas2)/2));
            rhos = nan(length(thetas),length(thetas2));
            
            for subch=1:length(oscChs{ch})
                for j = 1:length(thetas)
                    for o = 1:length(thetas2)
                        thetaRange1 = wrapTo2Pi(([thetas(j)+deg2rad(phaseBinWidth/2),thetas(j)-deg2rad(phaseBinWidth/2)]));
                        thetaRange2 = wrapTo2Pi(([thetas2(o)+deg2rad(phaseBinWidth/2),thetas2(o)-deg2rad(phaseBinWidth/2)]));

                        if thetaRange1(2) > thetaRange1(1) && thetaRange2(2) > thetaRange2(1)% Both Ranges crosses 0
                            thetaCount = length(find((data{ch,subch}<thetaRange1(1) | data{ch,subch}>thetaRange1(2)) & (data{end,1}<thetaRange2(1) | data{end,1}>thetaRange2(2))));
                        elseif thetaRange1(2) > thetaRange1(1) && thetaRange2(2) < thetaRange2(1)% Range 1 crosses 0
                            thetaCount = length(find((data{ch,subch}<thetaRange1(1) | data{ch,subch}>thetaRange1(2)) & (data{end,1}<thetaRange2(1) & data{end,1}>thetaRange2(2))));
                        elseif thetaRange1(2) < thetaRange1(1) && thetaRange2(2) > thetaRange2(1)% Range 2 crosses 0
                            thetaCount = length(find((data{ch,subch}<thetaRange1(1) & data{ch,subch}>thetaRange1(2)) & (data{end,1}<thetaRange2(1) | data{end,1}>thetaRange2(2))));
                        else %Neither Range crosses 0
                            thetaCount = length(find((data{ch,subch}<thetaRange1(1) & data{ch,subch}>thetaRange1(2)) & (data{end,1}<thetaRange2(1) & data{end,1}>thetaRange2(2))));
                        end
                        thetaPercent = thetaCount / sum(~isnan(data{end,1}));

                        if plotPercent == true
                            rhos(j,o) = mean([rhos(j,o),thetaPercent],'omitnan');
                        else
                            rhos(j,o) = sum([rhos(j,o),thetaCount],'omitnan');
                        end
                    end
                end
            end
            %Export data into output
            % output{n,1} = legendTitles{n};
            % output{n,2} = rmmissing(data)';
            % output{n,3} = [rad2deg(thetas(1:end-1))'];
            % output{n,4} = rhos;
            output{n,ch}=rhos;
            
            subplot(ceil(sqrt(length(oscChs)-1)),ceil((length(oscChs)-1) / ceil(sqrt(length(oscChs)-1))),ch)
            h(n,ch)=heatmap(rhos);
            h(n,ch).XDisplayLabels = string(ceil(wrapTo180(rad2deg(thetas2))));
            h(n,ch).YDisplayLabels = string(ceil(wrapTo180(rad2deg(thetas))));
            ylabel('PVxPFC→MD Phase')
            xlabel(oscChLabels{end})
            title(sprintf('%s',oscChLabels{ch}));
            set(gcf,'Position',[0 500 200 200])
            if max(rhos,[],'all') > MAX_COUNT
                MAX_COUNT = max(rhos,[],'all');
            end
        end
    end
end
for i = 1:size(h,1)
    for j = 1:size(h,2)
        set(h(i,j), 'ColorLimits', [0, MAX_COUNT*0.2]);
    end
end

%legend(legendTitles,'Location','southoutside')

% if size(output,1) == 2 && performWatson%if there are only two comparisons, perform a Watson's U2 test
%     [p,U2]=watsons_U2_approx_p(output{1,2},output{2,2});
%     fprintf('\nWatson U2 Test between %s and %s: p = %.4f, U2 = %.3f\n',output{1,1},output{2,1},p,U2);
%     annotation('textbox', [0.82, 0.15, 0.1, 0.1], 'String', sprintf("p=%.3f",p),'EdgeColor','none');
% end

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
        CHANNELNAMES{i} = 'wH-PVxMD Phases';
    elseif strcmp(oscChs{i},'xH')
        CHANNELNAMES{i} = 'xH-PVxMD Phases';
    elseif strcmp(oscChs{i},'PV')
        CHANNELNAMES{i} = 'PVxPV Phases';
    elseif strcmp(oscChs{i},'MD')
        CHANNELNAMES{i} = 'MDxMD Phases';
    elseif strcmp(oscChs{i},'wH-PD')
        CHANNELNAMES{i} = 'wH-PVxMD Phase Difference';
    elseif strcmp(oscChs{i},'xH-PD')
        CHANNELNAMES{i} = 'xH-PVxMD Phase Difference';
    elseif strcmp(oscChs{i},'wH-L')
        CHANNELNAMES{i} = 'wH-PVxMD (LH) Phases';
    elseif strcmp(oscChs{i},'xH-L')
        CHANNELNAMES{i} = 'xH-PVxMD (LH) Phases';
    elseif strcmp(oscChs{i},'wH-R')
        CHANNELNAMES{i} = 'wH-PVxMD (RH) Phases';
    elseif strcmp(oscChs{i},'xH-R')
        CHANNELNAMES{i} = 'xH-PVxMD (RH) Phases';
    end
end
if isscalar(CHANNELNAMES)
    CHTITLE = CHANNELNAMES{1};
else
    CHTITLE = 'Phases';
end

if plotThreshold==true
    THRESHTITLE = "Threshold: ON";
else
    THRESHTITLE = "Threshold: off";
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

function data = getData(currentCH,plotPhases,combIDX,combStateIDX)

%CHANNELS - MANUALLY DEFINED, POTENTIALL ALLOW FOR MORE FLEXIBILITY IN THE FUTURE

%From [[PVL,MDL,~,~,~];[PVL,PVR,~,~,~];[PVL,MDR,~,~,~];[PVR,MDL,~,~,~];[MDL,MDR,~,~,~];[PVR,MDR,~,~,~]]
wHLeft = 1;
PV = 2;
xHLeftPV = 3;
xHRightPV = 4;
MD = 5;
wHRight = 6;

%Select the appropriate rows from plotPhases
if strcmp(currentCH,'xH')
    data = [plotPhases(xHLeftPV,combIDX & combStateIDX) plotPhases(xHRightPV,combIDX & combStateIDX)];
elseif strcmp(currentCH,'wH')
    data = [plotPhases(wHLeft,(combIDX & combStateIDX)) plotPhases(wHRight,(combIDX & combStateIDX))];
elseif strcmp(currentCH,'PV')
    data = plotPhases(PV,combIDX & combStateIDX);
elseif strcmp(currentCH,'MD')
    data = plotPhases(MD,combIDX & combStateIDX);
elseif strcmp(currentCH,'xH-PD')
    data = wrapTo2Pi(abs(plotPhases(xHLeftPV,combIDX & combStateIDX)- plotPhases(xHRightPV,combIDX & combStateIDX)));
elseif strcmp(currentCH,'wH-PD')
    data = wrapTo2Pi(abs(plotPhases(wHLeft,combIDX & combStateIDX)- plotPhases(wHRight,combIDX & combStateIDX)));
elseif strcmp(currentCH,'xH-L')
    data = plotPhases(xHLeftPV,combIDX & combStateIDX);
elseif strcmp(currentCH,'wH-L')
    data = plotPhases(wHLeft,(combIDX & combStateIDX));
elseif strcmp(currentCH,'xH-R')
    data = plotPhases(xHRightPV,combIDX & combStateIDX);
elseif strcmp(currentCH,'wH-R')
    data = plotPhases(wHRight,(combIDX & combStateIDX));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

