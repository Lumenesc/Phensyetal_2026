%{

Rules for strategy identification:

C = Correct
E = Error
OC = Old Cue Selected (previously reinforced cue is part of choice)
NC = New Cue Selected (previously reinforced cue is NOT part of choice)

PER = Perseverative Block
REG = Regressive Block
NRB = New Behavior Block
REV = Reversal Block
CSS = Correct Shifted-Strategy Block

If first trial E + OC - immediately consider as PER Block

If following REG to C/E + NS -> NRB Block

if 4 OC in a row including 2+ errors -> REG Block
if 4 NC in a row including 2+ errors -> REV Block
if in 4 trials in a row both NC and OC including 2+ errors - NRB
if 3/4 Correct trials in a row -> CSS

Once in REG Block - do not break until NC trial or 3/4 correct
Once in REV Block - do not break until OC trial or 3/4 Correct
Once in a NRB block - do not break unless 2/4 Errors
Once in NSS - do not break unless 2/4 Errors

%}

%Global variables
C = 0;
E = 1;
NC = 0;
OC = 1;
PER = 0; REG = 1; NRB = 2; REV = 3; CSS = 4; NRB_OC = 5; NRB_NC = 6;
clear combStrat;
for k=1:size(timepoints,1)
    
    clear stratArr
    stratArr = zeros(size(timepoints{k,4},1),1);
    if timepoints{k,4} {1,7} == OC %OC Selected - should be default following IA
        stratArr(1) = PER;
    elseif timepoints{k,4} {1,7} == NC %New cue selected on first trial - Note: this should not happen but mice are unpredictable
        stratArr(1) = NRB;
    end
    
    for i=2:size(timepoints{k,4},1) %Move along mouse's trials in blocks of 4
        
        %if mean([timepoints{k,4}{i:end,6}]) <=0.2 %At anytime if the remaining trials are >80% correct
        
        %NOTE: First section rates last 8 trials as either CSS (Correct) or
        %NRB (Errors) specifically - necessary to avoid running out of blocks of 4
        if i>=(size(timepoints{k,4},1)-8) &&  timepoints{k,4}{i,6} == C %Last 8 correct trials are CSS
                stratArr(i) = CSS;
        elseif i>=(size(timepoints{k,4},1)-8) &&  timepoints{k,4}{i,6} == E %Last 8 correct trials are CSS
                stratArr(i) = NRB;
        else
        %{
        if stratArr(i-1) == CSS %Currently in CSS Block
            if sum([timepoints{k,4}{i:i+3,6}]) <=1 %Still in streak, continues CSS Block
                stratArr(i:i+3) = CSS;
            elseif sum([timepoints{k,4}{i:i+3,6}]) >=2 %2+ errors in a block of 4 breaks CSS Block
                
                if sum([timepoints{k,4}{i:i+3,7}]) == 4 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all OC and >2 errors
                     stratArr(i:i+3) = REG;
                elseif sum([timepoints{k,4}{i:i+3,7}]) == 0 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all NC and >2 errors
                     stratArr(i:i+3) = REV;
                elseif 0 > sum([timepoints{k,4}{i:i+3,7}]) < 4  %Next 4 trials are a mix of OC and NC
                     stratArr(i:i+3) = NRB;
                end
            end
         %}   
        if stratArr(i-1) == PER %Currently in PER Block
            
            
            if timepoints{k,4} {i,7} == OC %Old cue selected - continues PER Block
                stratArr(i:i+3) = PER;
            elseif timepoints{k,4} {i,7} == NC %New cue selected - breaks PER Block
                %if sum([timepoints{k,4}{i:i+3,6}]) <= 1 %If breaking block - 3/4 correct transitions NRB
                     %stratArr(i:i+3) = NRB;  
                if sum([timepoints{k,4}{i:i+3,7}]) >= 3 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all OC and >2 errors  TEMP SETTING TO REG, ORIGINALLY SET TO PER 7/10/23 & TEMP SET sum([timepoints{k,4}{i:i+3,6}]) >= 2 to sum([timepoints{k,4}{i:i+3,6}]) >= 1 07/23/23
                     stratArr(i:i+3) = PER;
                elseif sum([timepoints{k,4}{i:i+3,7}]) == 0 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all NC and >2 errors
                     stratArr(i:i+3) = REV;
                else %if 0 > sum([timepoints{k,4}{i:i+3,7}]) < 4  %Next 4 trials are a mix of OC and NC
                     stratArr(i:i+3) = NRB;
                end
            end
        elseif stratArr(i-1) == REG %Currently in REG Block
            
            if timepoints{k,4} {i,7} == OC %Old cue selected - continues REG Block
                stratArr(i:i+3) = REG;
            elseif timepoints{k,4} {i,7} == NC %New cue selected - breaks REG Block
                %if sum([timepoints{k,4}{i:i+3,6}]) <= 1 %If breaking block - 3/4 correct transitions NRB
                     %stratArr(i:i+3) = NRB;  
                if sum([timepoints{k,4}{i:i+3,7}]) >= 3 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all OC and >2 errors
                     stratArr(i:i+3) = REG;
                elseif sum([timepoints{k,4}{i:i+3,7}]) == 0 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all NC and >2 errors
                     stratArr(i:i+3) = REV;
                else %if 0 > sum([timepoints{k,4}{i:i+3,7}]) < 4  %Next 4 trials are a mix of OC and NC
                     stratArr(i:i+3) = NRB;
                end
            end
                
        elseif stratArr(i-1) == REV %Currently in REV Block
            
            if timepoints{k,4} {i,7} == NC %New cue selected - continues REV Block
                stratArr(i:i+3) = REV;
            elseif timepoints{k,4} {i,7} == OC %Old cue selected - breaks REV Block
                %if sum([timepoints{k,4}{i:i+3,6}]) <= 1 %If breaking block - 3/4 correct transitions NRB
                     %stratArr(i:i+3) = NRB;  
                if sum([timepoints{k,4}{i:i+3,7}]) == 4 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all OC and >2 errors
                     stratArr(i:i+3) = REG;
                elseif sum([timepoints{k,4}{i:i+3,7}]) <= 1 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all NC and >2 errors
                     stratArr(i:i+3) = REV;
                else %if 0 > sum([timepoints{k,4}{i:i+3,7}]) < 4 %Next 4 trials are a mix of OC and NC
                     stratArr(i:i+3) = NRB;
                end
            end
            
        elseif stratArr(i-1) == NRB || stratArr(i-1) == NRB_OC || stratArr(i-1) == NRB_NC %Currently in NRB Block
            
                if sum([timepoints{k,4}{i:i+3,7}]) == 4 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all OC and >2 errors
                     stratArr(i:i+3) = REG;
                elseif sum([timepoints{k,4}{i:i+3,7}]) == 0 && sum([timepoints{k,4}{i:i+3,6}]) >= 2 %Next 4 trials are all NC and >2 errors
                     stratArr(i:i+3) = REV;
                elseif 0 > sum([timepoints{k,4}{i:i+3,7}]) < 4 && timepoints{k,4}{i,7} == NC %Next 4 trials are a mix of OC and NC
                     stratArr(i:i+3) = NRB;
                elseif 0 > sum([timepoints{k,4}{i:i+3,7}]) < 4 && timepoints{k,4}{i,7} == OC %Next 4 trials are a mix of OC and NC
                     stratArr(i:i+3) = NRB;
                end
        end
        end
    end
    
    combStrat{k} = stratArr;
    
end

%Combine mice together with padding - aligned to end of behavior
clear temp
for i=1:size(timepoints,1)
temp{i,1}=combStrat{1,i};  %EDIT THIS
end
maxl = max(cellfun(@numel, {temp{:,1}}));
for i=1:size(timepoints,1)
temp{i,2} = padarray(temp{i,1},maxl-length(temp{i,1}),NaN,'pre');
end

combStrat{2,1} = [temp{:,2}];

%Combine mice together with padding - aligned to start of behavior
clear temp
for i=1:size(timepoints,1)
temp{i,1}=combStrat{1,i};  %EDIT THIS
end
maxl = max(cellfun(@numel, {temp{:,1}}));
for i=1:size(timepoints,1)
temp{i,2} = padarray(temp{i,1},maxl-length(temp{i,1}),NaN,'post');
end

combStrat{2,2} = [temp{:,2}];

%Concatenate all trials
clear temp Sidx
temp = double.empty;
for i=1:size(timepoints,1)
temp = [temp;combStrat{1,i}];  %EDIT THIS
end

combStrat{2,3} = temp;
Sidx = temp;
%{
%Create index of C/E, p.C/p.E, Mouse ID, gender, etc
clear Aidx
n=1;
for k=1:length(timepoints)
    for i=1:size(timepoints{k,3},1)
        
        Aidx{1,1} (n,1) = k;
        Aidx{1,2} (n,1) = timepoints{k,3}{i,6};
        Aidx{1,6} (n,1) = i;
        Aidx{1,7} (n,1) = size(timepoints{k,3},1);
        Aidx{1,8} (n,1) = 0; %0=IA  1=RS
        Aidx{1,10} (n,1) = timepoints{k,7}(1);      %Rewarded Cue 1=S, 2=L, 3=G, 4=C
        
        if timepoints{k,2} == 'F'
            Aidx{1,5} (n,1)= 1;
        elseif timepoints{k,2} == 'M'
            Aidx{1,5} (n,1) = 0;
        end
        
        %C/E
        if i == 1
            Aidx{1,3} (n,1) = timepoints{k,3}{i,6};
        else
            Aidx{1,3} (n,1) = timepoints{k,3}{i-1,6};
        end
        
        if size(timepoints{k,3},1) < 11
            Aidx{1,4} (n,1) = 0;
        elseif size(timepoints{k,3},1) < 21
            Aidx{1,4} (n,1) = 1;
        else
            Aidx{1,4} (n,1) = 2;
        end
        
        
        if i <= 5 %F5 of IA
            Aidx{1,9} (n,1) = 1;
        elseif  i >= size(timepoints{k,3},1)-10 %IA Strk
            Aidx{1,9} (n,1) = 3;
        else %Intermediate IA
            Aidx{1,9} (n,1) = 2;
        end
        
        n = n+1;
    end
end
for k=1:length(timepoints)
    for i=1:size(timepoints{k,4},1)
        
        Aidx{1,1} (n,1) = k;
        Aidx{1,2} (n,1) = timepoints{k,4}{i,6};
        Aidx{1,6} (n,1) = i;
        Aidx{1,7} (n,1) = size(timepoints{k,4},1);
        Aidx{1,8} (n,1) = 1; %0=IA  1=RS
        Aidx{1,10} (n,1) = timepoints{k,7}(2);      %Rewarded Cue 1=S, 2=L, 3=G, 4=C
        
        if timepoints{k,2} == 'F'
            Aidx{1,5} (n,1)= 1;
        elseif timepoints{k,2} == 'M'
            Aidx{1,5} (n,1) = 0;
        end
            
        %C/E
        if i == 1
            Aidx{1,3} (n,1) = timepoints{k,4}{i,6};
        else
            Aidx{1,3} (n,1) = timepoints{k,4}{i-1,6};
        end
        
        if size(timepoints{k,4},1) < 15
            Aidx{1,4} (n,1) = 0;
        elseif size(timepoints{k,4},1) < 30
            Aidx{1,4} (n,1) = 1;
        else
            Aidx{1,4} (n,1) = 2;
        end
        
        if i <= 5 %F5 of RS
            Aidx{1,9} (n,1) = 1;
        elseif i >= size(timepoints{k,4},1)-10 %RS Strk
            Aidx{1,9} (n,1) = 3;
        else %Intermediate RS
            Aidx{1,9} (n,1) = 2;
        end
        
        n = n+1;
    end
end
    
%}

clear C E NC OC REG NRB REV PER CSS stratArr i k temp maxl n

