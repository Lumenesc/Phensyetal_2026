
function [sigMatrixCmplx,metadata] = generateSigMatrixComplex(runName,sigMatrix,phaseoffsetsIndex,horzSmoothWindow,vertSmoothWindow,peakThreshold)

%% Function Description
%{

generateSigMatrixComplex v1 - AJP 11/21/2024

This function takes a sigMatrix with phase-shifts and converts it to a
matrix containing a complex vector consisting of the peak correlation and
it's phase for each timepoint across the recording session.

runName - User supplied run name

sigMatrix - generated from plot_temposcript_XCORR utilizing
temposcript_XCORR_PS

phaseoffsetsIndex - this is the specific phase angle for each row in the
sigMatrix{k,3}, this is defined by totalSS and xCorrShft in 
temposcript_XCORR_PS and should be formatted as such:
[-totalSS:xCorrShft:totalSS]. eg: [-360:22.5:360].

horzSmoothWindow - optional parameter used to smooth the correlation values
across time, this can help reduce the impact of artifacts, but may lead to
the loss of brief periods of synchrony

vertSmoothWindow - optional parameter used to smooth the correlation values
across phase shifts, this may help identify peak correlations but there is
currently not a strong motivation to use this parameter

peakThreshold - optional parameter. This sets a requirement that a peak
correlation is only considered if it is above this threshold (defined from
0 to 100) which is applied as a percentile of all correlation values across
phase shifts within each timepoint. If left blank, defaults to 0.

output:
    sigMatrixCmplx is the converted form of sigMatrix, note that now
    sigMatrix{k,3} is one row instead of many and is a complex array
    
    metadata is the parameter info used to perform conversion



%}
%% Set up variables and metadata



fprintf('Generating sigMatrixCmplx for %s...',runName);

if size(sigMatrix{1,3},1) ~= length(phaseoffsetsIndex)
    error('Phase Offsets Indices do not match the number of phase-offsets in sigMatrix');
else
    MF_radRange=(deg2rad(phaseoffsetsIndex));
end

if isempty(peakThreshold)
    peakThreshold = 0; %No Threshold for defining a peak phase-preference
end

metadata = struct('sigMatrixMeta',sigMatrix{1,2},'NumberofPhaseOffsets',size(sigMatrix{1,3},1),'PhaseOffsetsIndex_Degrees',phaseoffsetsIndex,'PeakPhasePreferenceThreshold',peakThreshold,'horzSmoothWindow',horzSmoothWindow,'vertSmoothWindow',vertSmoothWindow);

%% Apply smoothing algorithms
%Smoothing Algorithm - Uses Gaussian Window

%Horizontal smoothing will smooth across timepoints - useful for reducing
%the impact of noisy artifacts
if ~isempty(horzSmoothWindow) && horzSmoothWindow ~= 0
    for k=1:size(sigMatrix,1)
        sigMatrix{k,3} = smoothdata(sigMatrix{k,3},2,'gaussian',horzSmoothWindow); %Applies a gaussian smoothing function within each row
        sigMatrix{k,4} = smoothdata(sigMatrix{k,4},2,'gaussian',horzSmoothWindow); %Applies a gaussian smoothing function within each row - Control Comparison
    end
end

%Vertical smoothing will smooth across the phase offsets - in practice I
%have not found this particularily useful
if ~isempty(vertSmoothWindow) && vertSmoothWindow ~= 0
    for k=1:size(sigMatrix,1)
        sigMatrix{k,3} = smoothdata(sigMatrix{k,3},1,'gaussian',vertSmoothWindow); %Applies a gaussian smoothing function within each column
        sigMatrix{k,4} = smoothdata(sigMatrix{k,4},1,'gaussian',vertSmoothWindow); %Applies a gaussian smoothing function within each column - Control Comparison
    end
end

%% Convert sigMatrix to a complex vector array and save to sigMatrixCmplx - main output


%Identify phase-preference per sample and then convert to radians and
%enable looping
%%{
clear sigMatrixCmplx
for k=1:size(sigMatrix,1)

    %Update console with iteration
    if k==1
        fprintf('%d/%d',k,size(sigMatrix,1));
    else
        fprintf(repmat('\b', 1, length(sprintf('%d/%d',k-1,size(sigMatrix,1)))));
        fprintf('%d/%d',k,size(sigMatrix,1));
    end

    %Build the sigMatrixCmplx
    sigMatrixCmplx{k,1}=sigMatrix{k,1}; %Animal ID is the same

    for tp=1:size(sigMatrix{k,3},2) %Iterate through each timepoint

        %Generate the complex vector for the signal channel comparisons
        clear pks idx idx2
        [pks,idx] = findpeaks(sigMatrix{k,3}(:,tp)); %Identify phase-offsets with a peak correlation
        [~,idx2] = max(pks); %Select the peak with the strongest correlation - this is the phase preference

        if isempty(idx) %If there is NO peak correlation (rare) generate a NaN value for that timepoint
            sigMatrixCmplx{k,2}(tp) = complex(NaN,NaN);
        elseif max(pks) < prctile(sigMatrix{k,3}(:,tp),peakThreshold)  %Peaks must be > user selected threshold more than the other phases - this NaNs out 'weak' phase preferences
            sigMatrixCmplx{k,2}(tp) = complex(NaN,NaN);
        else %There is a valid peak phase
            sigMatrixCmplx{k,2}(tp) = sigMatrix{k,3}(idx(idx2),tp) * complex(cos(MF_radRange(idx(idx2))),sin(MF_radRange(idx(idx2)))); %Transform r-value and phase and r-value into a complex number
        end

        %Perform same computation for the control channel comparison 
        clear pks idx idx2
        [pks,idx] = findpeaks(sigMatrix{k,4}(:,tp)); %Identify phase-offsets with a peak correlation
        [~,idx2] = max(pks); %Select the peak with the strongest correlation - this is the phase preference

        if isempty(idx) %If there is NO peak correlation (rare) generate a NaN value for that timepoint
            sigMatrixCmplx{k,3}(tp) = complex(NaN,NaN);
        elseif max(pks) < prctile(sigMatrix{k,4}(:,tp),peakThreshold)  %Peaks must be > user selected threshold more than the other phases - this NaNs out 'weak' phase preferences
            sigMatrixCmplx{k,3}(tp) = complex(NaN,NaN);
        else %There is a valid peak phase
            sigMatrixCmplx{k,3}(tp) = sigMatrix{k,4}(idx(idx2),tp) * complex(cos(MF_radRange(idx(idx2))),sin(MF_radRange(idx(idx2)))); %Transform r-value and phase and r-value into a complex number
        end


    end

end

fprintf(repmat('\b', 1, length(sprintf('%d/%d',k,size(sigMatrix,1)))))
fprintf('Complete \n');

end