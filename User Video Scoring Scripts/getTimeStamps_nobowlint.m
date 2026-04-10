function [stamps,approaches] = getTimeStamps_nobowlint(filename,framerate,frameonsets)

clear A B C
[~,~,A]  = xlsread(filename);
A = string(A);
B = erase(A,"[");
B = erase(B,"]");
B = cellfun(@str2num,B,'UniformOutput',false);
for i=1:size(B,1)
    C{i,1} = cell2mat(B(i,:));
end

clear stamps
stamps(:,1) = C{1,1}(1:2:end)';
stamps(:,2) = C{2,1}(1:end)';
if size(C,1) == 4 %Error and Correct Separated
    stamps(:,3) = sort([C{3,1}(1:end)';C{4,1}(1:end)'],'ascend');
else
    stamps(:,3) = C{3,1}(1:end)';
end
stamps(:,4) = C{1,1}(2:2:end)';
stamps(1:end-1,5) = C{1,1}(3:2:end)'-(1*framerate);
stamps(end,5) = stamps(end,4)+(15*framerate);

stamps = floor(frameonsets(stamps));
approaches = [];

end
