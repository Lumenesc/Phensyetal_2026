%Used to change timepoints sample rate from 1second (default) to matching data sample rate to analyze
%fs is samples per second

%AJPhensy 03/2023

function [timepoints] = convTimepointFS(tpIn,fs)

for k=1:size(tpIn,1)
    for stage=2:3
        for i=1:size(tpIn{k,stage},1)
            for j=1:5
                tpIn{k,stage}{i,j}= floor(tpIn{k,stage}{i,j} * fs);
            end
        end
    end
end


timepoints = tpIn;

end