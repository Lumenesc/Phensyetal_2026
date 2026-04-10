%Used when timepoints matrix is missing a final value for ITI End
%This is caused since ITI End is calculated from T.Start of next trial, but
%on final trial there is no next.

%Padding is decided by 'length' in seconds

function [timepoints] = padITIEnd(tpIn,length)

for k=1:size(tpIn,1)
    if isempty(tpIn{k,2}{end,5})
        tpIn{k,2}{end,5} = tpIn{k,2}{end,4}+length;
    end
    if isempty(tpIn{k,3}{end,5})
        tpIn{k,3}{end,5} = tpIn{k,3}{end,4}+length;
    end
end

timepoints = tpIn;

end