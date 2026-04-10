function [metaIA metaRS] = getMetadata(filenameIA,filenameRS)

clear A B C
[~,~,A]  = xlsread(filenameIA);
A = string(A);

[~,~,B]  = xlsread(filenameRS);
B = string(B);

clear metaIA metaRS
%Get C/E column
metaIA(:,1) = str2double(A(:,2));
%Convert L to 0 and R to 1
for i=1:size(A,1)
    if strcmp(A(i,3),'L')
        metaIA(i,3) = 0;
    elseif strcmp(A(i,3),'R')
        metaIA(i,3) = 1;
    else
        metaIA(i,3) = NaN;
    end
end

%RS Metadata
%Get C/E Column
metaRS(:,1) = str2double(B(:,2));

for i=1:size(B,1)
    if strcmp(B(i,3),'L')
        metaRS(i,3) = 0;
    elseif strcmp(B(i,3),'R')
        metaRS(i,3) = 1;
    else
        metaRS(i,3) = NaN;
    end
end

%Old Strat vs New Strat (RS Only)
%Pair1: SC/LG, Pair2: SG/LC
for i=1:size(B,1)
    if strcmp(B(i,4),"1") && ...
         ((strcmp(A(1,5),"L") && strcmp(B(1,5),"G")) || ...
                (strcmp(A(1,5),"G") && strcmp(B(1,5),"L")) || ...
                (strcmp(A(1,5),"S") && strcmp(B(1,5),"C")) ||...
                (strcmp(A(1,5),"C") && strcmp(B(1,5),"S")))
        if strcmp(B(i,2),"0")
            metaRS(i,2) = 1;
        else 
            metaRS(i,2) = 0;
        end

    elseif strcmp(B(i,4),"2") && ...
         ((strcmp(A(1,5),"L") && strcmp(B(1,5),"C")) || ...
                (strcmp(A(1,5),"C") && strcmp(B(1,5),"L")) || ...
                (strcmp(A(1,5),"S") && strcmp(B(1,5),"G")) ||...
                (strcmp(A(1,5),"G") && strcmp(B(1,5),"S")))
        
        if strcmp(B(i,2),"0")
            metaRS(i,2) = 1;
        else 
            metaRS(i,2) = 0;
        end
    else
        if strcmp(B(i,2),"0")
            metaRS(i,2) = 0;
        else 
            metaRS(i,2) = 1;
        end
    end
end

%Stay/Go or P.Selected
%0: Selected same cue last presented (Stay), 1: Did not select same cue (Go)

for i=1:size(B,1)
    pair = str2double(B(i,4));
    if ~isempty(find(str2double(B(1:i-1,4))==pair,1,'last'))
        idx = find(str2double(B(1:i-1,4))==pair,1,'last');
        if B(i,2) == B(idx,2)
             metaRS(i,4) = 0;
        else
             metaRS(i,4) = 1;
        end
    else
        idx = find(str2double(A(1:end,4))==pair,1,'last');
        if metaRS(i,2) == str2double(A(idx,2))
             metaRS(i,4) = 1;
        else
             metaRS(i,4) = 0;
        end
    end
end

for i=5:size(B,1)
    for j=1:4
        if B(i,4) == B(i-5+j,4)
            if B(i,2) == B(i-5+j,2)
                metaRS(i,4) = 0;
            else
                metaRS(i,4) = 1;
            end
        end
    end
end
    

end