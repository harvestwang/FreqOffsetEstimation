function freqOffsetEst = crossCorrFreqEstimate(Seq1, Seq2, modRate,  interval)


% crossCorrSeq1 = decRxSigCorrect1(1 : syncLen);
% crossCorrSeq2 = decRxSigCorrect1(syncLen+dataLen+1 : frameLen);

    crossCorrRes = mean(Seq2 .* conj(Seq1));
    freqOffsetEst = angle(crossCorrRes)/(2*pi)*modRate/interval;

end

