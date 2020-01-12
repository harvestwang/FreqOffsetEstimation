function freqOffsetEst = newFreqEstimate(dpDephasePre, dpDephasePost, interval)
        % Fisrt frequency offset estimation and correction
        
%         k = (1 : syncLen-1)';
%         w1 = 6 .* k .* (syncLen-k)/(syncLen*(syncLen-1)); % window
%         
%         selfCorrPreSeq1 = dephasePre(1 : syncLen-alpha);
%         selfCorrPreSeq2 = dephasePre(alpha+1 : syncLen);
%         selfCorrPreRes = mean(w1 .* selfCorrPreSeq2 .* conj(selfCorrPreSeq1));
%         firFreqOffsetEstTemp(i, time) = angle(selfCorrPreRes)/(2*pi)*modRate/alpha;
        
%         dephaseRxCorrected = dephaseRx .* ...
%             exp(-1j*2*pi*(0:frameLen-1)'*firFreqOffsetEstTemp(i, time)/clkFreq);
%         rxSyncPreCorrect = dephasePre .* ...
%             exp(-1j*2*pi*(0:syncLen-1)'*firFreqOffsetEstTemp(i, time)/clkFreq);
%         rxSyncPostCorrect = dephasePost .* ...
%             exp(-1j*2*pi*(syncLen+dataLen:frameLen-1)'*firFreqOffsetEstTemp(i, time)/clkFreq);
        
        % Second frequency offset estimation and correction
%         crossCorrSeq1 = dephaseRxCorrected(1 : syncLen);
%         crossCorrSeq2 = dephaseRxCorrected(syncLen+dataLen+1 : frameLen);
%         interval = syncLen + dataLen;
        
%         secFreqOffsetEstTemp(i, time) = crossCorrFreqEstimate(crossCorrSeq1, ...
%             crossCorrSeq2, modRate, interval);
        
%         dephaseRxCorrected1 = dephaseRx .* ...
%             exp(-1j*2*pi*(0:frameLen-1)'*dftFreqOffsetEstTemp(i, time)/clkFreq);
        
%         crossCorrSeq1 = dephaseRxCorrected1(1 : syncLen);
%         crossCorrSeq2 = dephaseRxCorrected1(syncLen+dataLen+1 : frameLen);
%         interval = syncLen + dataLen;
        
%         secFreqOffsetEstTemp1(i, time) = crossCorrFreqEstimate(crossCorrSeq1, ...
%             crossCorrSeq2, modRate, interval);
        
%         newFreqOffsetEstTemp(i, time) = firFreqOffsetEstTemp(i, time) + secFreqOffsetEstTemp(i, time);
%         newFreqOffsetEstTemp1(i, time) = dftFreqOffsetEstTemp(i, time) + secFreqOffsetEstTemp1(i, time);

        % second frequency offset correction
        %         decRxSigCorrect2 = decRxSigCorrect1 .* ...
        %             exp(-1j*2*pi*(0:syncLen-1)'*secFreqOffsetEstTemp(i, time)/clkFreq);
end

