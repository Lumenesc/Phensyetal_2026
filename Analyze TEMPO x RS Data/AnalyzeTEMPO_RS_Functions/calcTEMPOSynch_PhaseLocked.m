
function synchronyValues = calcTEMPOSynch_PhaseLocked(sigMatrixCmplx,timepoints,phaserange,trSecTimes,normON)

clear zstgValues zstgRngValues zstgValuesSmooth

[zstgValues(:,1),zstgRngValues{1,1},zstgValuesSmooth(:,1),~] = get_stgPhaseRangeSynchbyTrial_norm(sigMatrixCmplx,timepoints,2,phaserange,trSecTimes,normON); 
[zstgValues(:,2),zstgRngValues{1,2},zstgValuesSmooth(:,2),stgRngTimes] = get_stgPhaseRangeSynchbyTrial_norm(sigMatrixCmplx,timepoints,3,phaserange,trSecTimes,normON);

synchronyValues = {zstgValues,zstgRngValues,zstgValuesSmooth,stgRngTimes};

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stgValues,stgRngValues,stgValuesSmooth,stgRngTimes] = get_stgPhaseRangeSynchbyTrial_norm(sigSHFCompPeaks,timepoints,taskstg,phrangedeg,trSecTimes,normON) 
PREDEC = 1; DIG = 2; OUTCOME = 3; POSTOUT = 4; ITIST = 5; ITIEND = 6; MAXPREDEC = 7;
STAMP1 = 1; OFFSET1 = 2; STAMP2 = 3; OFFSET2 = 4;

    %stgRngValues = NaN; %When Disabled to allow variable ranges of stages
    
    %phrange is expected in degrees

    fprintf("Function Called: get_stgPhaseRangeSynchbyTrial\n");

    %Generate an IDX_pair for the next-pair presentation task-phases
    if taskstg==3
    for k=1:size(timepoints,1)
        for tr=1:size(timepoints{k,taskstg},1)
            if  timepoints{k,taskstg}{tr,6}==0 && timepoints{k,taskstg}{tr,7}==1 %Correct and OS (ReinfPair)
                IDX_pair{k,1}(tr,1) = 1;
            elseif timepoints{k,taskstg}{tr,6}==1 && timepoints{k,taskstg}{tr,7}==1 %Error and OS (NewPair)
                IDX_pair{k,1}(tr,1) = 2;
            elseif timepoints{k,taskstg}{tr,6}==0 && timepoints{k,taskstg}{tr,7}==0 %Correct and NS (NewPair)
                IDX_pair{k,1}(tr,1) = 2;
            elseif timepoints{k,taskstg}{tr,6}==1 && timepoints{k,taskstg}{tr,7}==0 %Error and NS (ReinfPair)
                IDX_pair{k,1}(tr,1) = 1;
            end
        end
    end
    end

    for k=1:size(timepoints,1)
        fprintf("\nAnimal: %s. Trial: ",timepoints{k,1});
        for tr=1:size(timepoints{k,taskstg},1)
            stage = 0;
            fprintf("%d ",tr);
            

            
            %Previous Non-Conflict Post-Outcome 
            stage = stage+1;
            if tr==1 || taskstg==2 || sum(IDX_pair{k,1}(1:tr-1,1)==1)==0
                stgValues{k,1}(tr,stage) = NaN;
                stgRngValues{k,stage}(tr,:) =  nan(1,length(4*-10:4*10));
            else
                ptr = find(IDX_pair{k,1}(1:tr-1,1)==1,1,'last');
                stamp1 = trSecTimes{2,POSTOUT}(STAMP1); stamp2 = trSecTimes{2,POSTOUT}(STAMP2);
                preOffset = 4 * trSecTimes{2,POSTOUT}(OFFSET1); postOffset = 4 * trSecTimes{2,POSTOUT}(OFFSET2);
                tp1 = timepoints{k,taskstg}{ptr,stamp1}+preOffset;
                tp2 = timepoints{k,taskstg}{ptr,stamp2}+postOffset;
                stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
                stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

                tp1 = timepoints{k,taskstg}{ptr,3}-(4*10);
                tp2 = timepoints{k,taskstg}{ptr,3}+(4*10);
                stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);
            end
            
            %Previous Conflict Post-Outcome 
            stage = stage+1;
            if tr==1 || taskstg==2 || sum(IDX_pair{k,1}(1:tr-1,1)==2)==0
                stgValues{k,1}(tr,stage) = NaN;
                stgRngValues{k,stage}(tr,:) =  nan(1,length(4*-10:4*10));
            else
                ptr = find(IDX_pair{k,1}(1:tr-1,1)==2,1,'last');
                stamp1 = trSecTimes{2,POSTOUT}(STAMP1); stamp2 = trSecTimes{2,POSTOUT}(STAMP2);
                preOffset = 4 * trSecTimes{2,POSTOUT}(OFFSET1); postOffset = 4 * trSecTimes{2,POSTOUT}(OFFSET2);
                tp1 = timepoints{k,taskstg}{ptr,stamp1}+preOffset;
                tp2 = timepoints{k,taskstg}{ptr,stamp2}+postOffset;
                stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
                stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

                tp1 = timepoints{k,taskstg}{ptr,3}-(4*10);
                tp2 = timepoints{k,taskstg}{ptr,3}+(4*10);
                stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);
            end
         
            %Pre-Decision
            stage = stage+1;
            if timepoints{k,taskstg}{tr,2}-timepoints{k,taskstg}{tr,1} > 4*trSecTimes{2,MAXPREDEC}%CHANGED TO IF GREATER THAN 10SECONDS
                stamp1 = 2; stamp2 = trSecTimes{2,PREDEC}(STAMP2);
                preOffset = 4 * -trSecTimes{2,MAXPREDEC}; postOffset = 4 * trSecTimes{2,PREDEC}(OFFSET2);
                tp1 = timepoints{k,taskstg}{tr,stamp1}+preOffset;
                tp2 = timepoints{k,taskstg}{tr,stamp2}+postOffset;
            else
                stamp1 = trSecTimes{2,PREDEC}(STAMP1); stamp2 = trSecTimes{2,PREDEC}(STAMP2);
                preOffset = 4 * trSecTimes{2,PREDEC}(OFFSET1); postOffset = 4 * trSecTimes{2,PREDEC}(OFFSET2);
                tp1 = timepoints{k,taskstg}{tr,stamp1}+preOffset;
                tp2 = timepoints{k,taskstg}{tr,stamp2}+postOffset;
            end
            stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
            stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

            tp1 = timepoints{k,taskstg}{tr,2}-(4*10);
            tp2 = timepoints{k,taskstg}{tr,2}+(4*10);
            stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);

            %Dig
            stage = stage+1;
            stamp1 = trSecTimes{2,DIG}(STAMP1); stamp2 = trSecTimes{2,DIG}(STAMP2);
            preOffset = 4 * trSecTimes{2,DIG}(OFFSET1); postOffset = 4 * trSecTimes{2,DIG}(OFFSET2);
            tp1 = timepoints{k,taskstg}{tr,stamp1}+preOffset;
            tp2 = timepoints{k,taskstg}{tr,stamp2}+postOffset;
            stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
            stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

            tp1 = timepoints{k,taskstg}{tr,2}-(4*10);
            tp2 = timepoints{k,taskstg}{tr,2}+(4*10);
            stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);

            %Outcome
            stage = stage+1;
            stamp1 = trSecTimes{2,OUTCOME}(STAMP1); stamp2 = trSecTimes{2,OUTCOME}(STAMP2);
            preOffset = 4 * trSecTimes{2,OUTCOME}(OFFSET1); postOffset = 4 * trSecTimes{2,OUTCOME}(OFFSET2);
            tp1 = timepoints{k,taskstg}{tr,stamp1}+preOffset;
            tp2 = timepoints{k,taskstg}{tr,stamp2}+postOffset;
            stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
            stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];
            stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);

            %Post-Outcome
            stage = stage+1;
            stamp1 = trSecTimes{2,POSTOUT}(STAMP1); stamp2 = trSecTimes{2,POSTOUT}(STAMP2);
            preOffset = 4 * trSecTimes{2,POSTOUT}(OFFSET1); postOffset = 4 * trSecTimes{2,POSTOUT}(OFFSET2);
            tp1 = timepoints{k,taskstg}{tr,stamp1}+preOffset;
            tp2 = timepoints{k,taskstg}{tr,stamp2}+postOffset;
            stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
            stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

            tp1 = timepoints{k,taskstg}{tr,3}-(4*10);
            tp2 = timepoints{k,taskstg}{tr,3}+(4*10);
            stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);

            
            %ITI-First Half
            stage = stage+1;
            if tr<size(timepoints{k,taskstg},1)%Not last trial, there exists a next trial
            stamp1 = trSecTimes{2,ITIST}(STAMP1); stamp2 = trSecTimes{2,ITIST}(STAMP2);
            preOffset = 4 * trSecTimes{2,ITIST}(OFFSET1); postOffset = 4 * trSecTimes{2,ITIST}(OFFSET2);
            tp1 = timepoints{k,taskstg}{tr,stamp1}+preOffset;
            tp2 = timepoints{k,taskstg}{tr,stamp2}+postOffset;
            stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
            stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

            tp1 = timepoints{k,taskstg}{tr,4}-(4*0);
            tp2 = timepoints{k,taskstg}{tr,4}+(4*15);
            stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);
            else
                stgValues{k,1}(tr,stage) = NaN;
                stgRngValues{k,stage}(tr,:) = NaN;
            end

            %ITI-Second Half
            stage = stage+1;
            if tr<size(timepoints{k,taskstg},1)%Not last trial, there exists a next trial
            stamp1 = trSecTimes{2,ITIEND}(STAMP1); stamp2 = trSecTimes{2,ITIEND}(STAMP2);
            preOffset = 4 * trSecTimes{2,ITIEND}(OFFSET1); postOffset = 4 * trSecTimes{2,ITIEND}(OFFSET2);
            tp1 = timepoints{k,taskstg}{tr,stamp1}+preOffset;
            tp2 = timepoints{k,taskstg}{tr,stamp2}+postOffset;
            stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
            stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

            tp1 = timepoints{k,taskstg}{tr,5}-(4*15);
            tp2 = timepoints{k,taskstg}{tr,5}+(4*0);
            stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);
            stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);
            else
                stgValues{k,1}(tr,stage) = NaN;
                stgRngValues{k,stage}(tr,:) = NaN;
            end

            %Next-Trial Pre-Decision
            stage = stage+1;
            if tr<size(timepoints{k,taskstg},1) %Not last trial, there exists a next trial
                if timepoints{k,taskstg}{tr+1,2}-timepoints{k,taskstg}{tr+1,1} > 4*trSecTimes{2,MAXPREDEC}%CHANGED TO IF GREATER THAN 10SECONDS
                    stamp1 = 2; stamp2 = trSecTimes{2,PREDEC}(STAMP2);
                    preOffset = 4 * -trSecTimes{2,MAXPREDEC}; postOffset = 4 * trSecTimes{2,PREDEC}(OFFSET2);
                    tp1 = timepoints{k,taskstg}{tr+1,stamp1}+preOffset;
                    tp2 = timepoints{k,taskstg}{tr+1,stamp2}+postOffset;
                else
                    stamp1 = trSecTimes{2,PREDEC}(STAMP1); stamp2 = trSecTimes{2,PREDEC}(STAMP2);
                    preOffset = 4 * trSecTimes{2,PREDEC}(OFFSET1); postOffset = 4 * trSecTimes{2,PREDEC}(OFFSET2);
                    tp1 = timepoints{k,taskstg}{tr+1,stamp1}+preOffset;
                    tp2 = timepoints{k,taskstg}{tr+1,stamp2}+postOffset;
                end
                stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
                stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

                tp1 = timepoints{k,taskstg}{tr,2}-(4*10);
                tp2 = timepoints{k,taskstg}{tr,2}+(4*10);
                stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);

            else
                stgValues{k,1}(tr,stage) = NaN;
                stgRngValues{k,stage}(tr,:) = NaN;
            end
            
            
            %Next Non-Conflict Pre-Decision 
            stage = stage+1;
            PAIR = 1; %Non-Conflict
            if tr==size(timepoints{k,taskstg},1) || taskstg==2 || sum(IDX_pair{k,1}(tr+1:end,1)==PAIR)==0
                stgValues{k,1}(tr,stage) = NaN;
                stgRngValues{k,stage}(tr,:) =  nan(1,length(4*-10:4*10));
            else
                idx = find(IDX_pair{k,1}==PAIR); %Grab an index of all RS trials of PAIR
                ptr = idx(sum(IDX_pair{k,1}(1:tr,1)==PAIR)+1); %Identify which is the next PAIR trial
                if timepoints{k,taskstg}{ptr,2}-timepoints{k,taskstg}{ptr,1} > 4*trSecTimes{2,MAXPREDEC}%CHANGED TO IF GREATER THAN 10SECONDS
                    stamp1 = 2; stamp2 = trSecTimes{2,PREDEC}(STAMP2);
                    preOffset = 4 * -trSecTimes{2,MAXPREDEC}; postOffset = 4 * trSecTimes{2,PREDEC}(OFFSET2);
                    tp1 = timepoints{k,taskstg}{ptr,stamp1}+preOffset;
                    tp2 = timepoints{k,taskstg}{ptr,stamp2}+postOffset;
                else
                    stamp1 = trSecTimes{2,PREDEC}(STAMP1); stamp2 = trSecTimes{2,PREDEC}(STAMP2);
                    preOffset = 4 * trSecTimes{2,PREDEC}(OFFSET1); postOffset = 4 * trSecTimes{2,PREDEC}(OFFSET2);
                    tp1 = timepoints{k,taskstg}{ptr,stamp1}+preOffset;
                    tp2 = timepoints{k,taskstg}{ptr,stamp2}+postOffset;
                end
                stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
                stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

                tp1 = timepoints{k,taskstg}{ptr,3}-(4*10);
                tp2 = timepoints{k,taskstg}{ptr,3}+(4*10);
                stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);
            end
            
            %Next Conflict Pre-Decision 
            stage = stage+1;
            PAIR = 2; %Conflict
            if tr==size(timepoints{k,taskstg},1) || taskstg==2 || sum(IDX_pair{k,1}(tr+1:end,1)==PAIR)==0
                stgValues{k,1}(tr,stage) = NaN;
                stgRngValues{k,stage}(tr,:) =  nan(1,length(4*-10:4*10));
            else
                idx = find(IDX_pair{k,1}==PAIR); %Grab an index of all RS trials of PAIR
                ptr = idx(sum(IDX_pair{k,1}(1:tr,1)==PAIR)+1); %Identify which is the next PAIR trial
                if timepoints{k,taskstg}{ptr,2}-timepoints{k,taskstg}{ptr,1} > 4*trSecTimes{2,MAXPREDEC}%CHANGED TO IF GREATER THAN 10SECONDS
                    stamp1 = 2; stamp2 = trSecTimes{2,PREDEC}(STAMP2);
                    preOffset = 4 * -trSecTimes{2,MAXPREDEC}; postOffset = 4 * trSecTimes{2,PREDEC}(OFFSET2);
                    tp1 = timepoints{k,taskstg}{ptr,stamp1}+preOffset;
                    tp2 = timepoints{k,taskstg}{ptr,stamp2}+postOffset;
                else
                    stamp1 = trSecTimes{2,PREDEC}(STAMP1); stamp2 = trSecTimes{2,PREDEC}(STAMP2);
                    preOffset = 4 * trSecTimes{2,PREDEC}(OFFSET1); postOffset = 4 * trSecTimes{2,PREDEC}(OFFSET2);
                    tp1 = timepoints{k,taskstg}{ptr,stamp1}+preOffset;
                    tp2 = timepoints{k,taskstg}{ptr,stamp2}+postOffset;
                end
                stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
                stgRngTimes(stage,1:4) = [stamp1 stamp2 preOffset postOffset];

                tp1 = timepoints{k,taskstg}{ptr,3}-(4*10);
                tp2 = timepoints{k,taskstg}{ptr,3}+(4*10);
                stgRngValues{k,stage}(tr,:) = sigSHFCompPeaks{k,2}(:,tp1:tp2);
            end
            

            %Full Trial (Non-ITI)
            stage = stage+1;
            tp1 = timepoints{k,taskstg}{tr,1};
            tp2 = timepoints{k,taskstg}{tr,4};
            stgValues{k,1}(tr,stage) = get_phasesynch_norm(sigSHFCompPeaks{k,2},tp1,tp2,phrangedeg,normON);
            
        end

        %Incorporating a smoothing algorithm across trials
        stgValuesSmooth{k,1} = smoothdata(stgValues{k,1},1,'gaussian',10);

    end
    fprintf("\n");
end

function [synchOUT] = get_phasesynch_norm(cmplxdata,tp1,tp2,phrange,normON)

    %Find the timepoints with a phase within phrange and average synchrony

    baseline = [240 6000]; %Range used to calculate z-score

    window = tp2-tp1;
    
    phrange = wrapTo2Pi(deg2rad(phrange));
    
    
    %If only 2 range variables (single range)

    if length(phrange) < 4

        if phrange(1) > phrange(2) %Range crosses 0 degrees

            %First get avereage and std of full trial based on window size
            n=0;
            clear syncbywindow
            for i=baseline(1):window:baseline(2)%size(cmplxdata,2)-window
                n=n+1;
                cmplxseg = cmplxdata(i:i+window);
                syncbywindow(n) = mean(abs(cmplxseg(wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(1) | wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(2))),'omitnan');
            end
            %syncbywindow(isnan(syncbywindow))=0;
            zmean = mean(syncbywindow,"omitnan");
            zstd = std(syncbywindow,"omitnan");

            %Calculate z-scored peak synchrony
            cmplxseg = cmplxdata(tp1:tp2);
            synch = mean(abs(cmplxseg(wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(1) | wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(2))),'omitnan');
            if isnan(synch)
                zsynch = (0 - zmean) / zstd;
            else
                zsynch = (synch - zmean) / zstd;
            end

        else
        
            %First get avereage and std of full trial based on window size
            n=0;
            clear syncbywindow
            for i=baseline(1):window:baseline(2)%size(cmplxdata,2)-window
                n=n+1;
                cmplxseg = cmplxdata(i:i+window);
                syncbywindow(n) = mean(abs(cmplxseg(wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(1) & wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(2))),'omitnan');
            end
            syncbywindow(isnan(syncbywindow))=0;
            zmean = mean(syncbywindow,"omitnan");
            zstd = std(syncbywindow,"omitnan");

            %Calculate z-scored peak synchrony
            cmplxseg = cmplxdata(tp1:tp2);
            synch = mean(abs(cmplxseg(wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(1) & wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(2))),'omitnan');
            if isnan(synch)
                zsynch = (0 - zmean) / zstd;
            else
                zsynch = (synch - zmean) / zstd;
            end
        end
    
    %Otherwise there's 2 sets of ranges
    else 
        if phrange(1) > phrange(2) && phrange(3) > phrange(4) %Both ranges cross 0 degrees

            %First get avereage and std of full trial based on window size
            n=0;
            clear syncbywindow
            for i=baseline(1):window:baseline(2)
                n=n+1;
                cmplxseg = cmplxdata(i:i+window);
                syncbywindow(n) = mean(abs(cmplxseg((wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(1) | wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(2)) | ...
                                                (wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(3) | wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(4)) )),'omitnan');
            end
            syncbywindow(isnan(syncbywindow))=0;
            zmean = mean(syncbywindow,"omitnan");
            zstd = std(syncbywindow,"omitnan");

            %Calculate z-scored peak synchrony
            cmplxseg = cmplxdata(tp1:tp2);
            synch = mean(abs(cmplxseg((wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(1) | wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(2)) | ...
                                            (wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(3) | wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(4)) )),'omitnan');
            if isnan(synch)
                zsynch = (0 - zmean) / zstd;
            else
                zsynch = (synch - zmean) / zstd;
            end
        elseif phrange(1) > phrange(2) %First range crosses 0 degrees
            %First get avereage and std of full trial based on window size
            n=0;
            clear syncbywindow
            for i=baseline(1):window:baseline(2)
                n=n+1;
                cmplxseg = cmplxdata(i:i+window);
                syncbywindow(n) = mean(abs(cmplxseg((wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(1) | wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(2)) | ...
                                                (wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(3) & wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(4)) )),'omitnan');
            end
            syncbywindow(isnan(syncbywindow))=0;
            zmean = mean(syncbywindow,"omitnan");
            zstd = std(syncbywindow,"omitnan");

            %Calculate z-scored peak synchrony
            cmplxseg = cmplxdata(tp1:tp2);
            synch = mean(abs(cmplxseg((wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(1) | wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(2)) | ...
                                            (wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(3) & wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(4)) )),'omitnan');
            if isnan(synch)
                zsynch = (0 - zmean) / zstd;
            else
                zsynch = (synch - zmean) / zstd;
            end
        elseif phrange(3) > phrange(4) %Second range crosses 0 degrees
            %First get avereage and std of full trial based on window size
            n=0;
            clear syncbywindow
            for i=baseline(1):window:baseline(2)
                n=n+1;
                cmplxseg = cmplxdata(i:i+window);
                syncbywindow(n) = mean(abs(cmplxseg((wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(1) & wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(2)) | ...
                                                (wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(3) | wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(4)) )),'omitnan');
            end
            syncbywindow(isnan(syncbywindow))=0;
            zmean = mean(syncbywindow,"omitnan");
            zstd = std(syncbywindow,"omitnan");

            %Calculate z-scored peak synchrony
            cmplxseg = cmplxdata(tp1:tp2);
            synch = mean(abs(cmplxseg((wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(1) & wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(2)) | ...
                                            (wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(3) | wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(4)) )),'omitnan');
            if isnan(synch)
                zsynch = (0 - zmean) / zstd;
            else
                zsynch = (synch - zmean) / zstd;
            end
        else %Neither range crosses 0 degrees
            %First get avereage and std of full trial based on window size
            n=0;
            clear syncbywindow
            for i=baseline(1):window:baseline(2)
                n=n+1;
                cmplxseg = cmplxdata(i:i+window);
                syncbywindow(n) = mean(abs(cmplxseg((wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(1) & wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(2)) | ...
                                                (wrapTo2Pi(angle(cmplxdata(i:i+window))) > phrange(3) & wrapTo2Pi(angle(cmplxdata(i:i+window))) < phrange(4)) )),'omitnan');
            end
            syncbywindow(isnan(syncbywindow))=0;
            zmean = mean(syncbywindow,"omitnan");
            zstd = std(syncbywindow,"omitnan");

            %Calculate z-scored peak synchrony
            cmplxseg = cmplxdata(tp1:tp2);
            synch = mean(abs(cmplxseg((wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(1) & wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(2)) | ...
                                            (wrapTo2Pi(angle(cmplxdata(tp1:tp2))) > phrange(3) & wrapTo2Pi(angle(cmplxdata(tp1:tp2))) < phrange(4)) )),'omitnan');
            if isnan(synch)
                zsynch = (0 - zmean) / zstd;
            else
                zsynch = (synch - zmean) / zstd;
            end
        end
    end

    if normON==true
        synchOUT = zsynch;
    else
        if isnan(synch)
            synchOUT = 0;
        else
            synchOUT = synch;
        end
    end


end