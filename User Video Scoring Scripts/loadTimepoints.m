
k = size(timepoints,1)+1;
tp_row = cell.empty(0,3);

tp_row{1,1} = 'H2C3_5';


%Load tank file in order to correct for dropped frames and convert to
%seconds
if exist('tank_path','var')
    [tank_path] = uigetdir(tank_path,'Synapse Tank Folder');
else
    [tank_path] = uigetdir('Title','Synapse Tank Folder');
end
[data] =  TDTbin2mat(tank_path);
frameonsets = data.epocs.Cam1.onset;

%First Load IA
if exist('path','var')
    [file,path] = uigetfile({'*.xlsx';'*.xls';'*.xlsm'}, ...
        'IA VideoTracker Data',path);
else
    [file,path] = uigetfile({'*.xlsx';'*.xls';'*.xlsm'}, ...
        'IA VideoTracker Data');
end
[stamps,~] = getTimeStamps_nobowlint([path,file],20,frameonsets);
tp_row{1,2} = num2cell(stamps);

%Then Load RS
[file,path] = uigetfile({'*.xlsx';'*.xls';'*.xlsm'}, ...
'RS VideoTracker Data',path);
[stamps,~] = getTimeStamps_nobowlint([path,file],20,frameonsets);
tp_row{1,3} = num2cell(stamps);

clear metaIAfile metaRSfile metaIA metaRS
[file,path] = uigetfile({'*.xlsx';'*.xls';'*.xlsm'}, ...
'IA MetaData',path);
metaIAfile = [path,file];
[file,path] = uigetfile({'*.xlsx';'*.xls';'*.xlsm'}, ...
'RS MetaData',path);
metaRSfile = [path,file];

[metaIA,metaRS] = getMetadata(metaIAfile,metaRSfile);
tp_row{1,2}(:,6:8) = num2cell(metaIA);
tp_row{1,3}(:,6:9) = num2cell(metaRS);

timepoints(k,:) = tp_row;

clear stamps metaRS metaIA metaIAfile metaRSfile k file