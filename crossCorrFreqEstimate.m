function freqOffsetEst = crossCorrFreqEstimate( crossCorrSeq1, crossCorrSeq2, modRate,  interval)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% crossCorrSeq1 = decRxSigCorrect1(1 : syncLen);
% crossCorrSeq2 = decRxSigCorrect1(syncLen+dataLen+1 : frameLen);

crossCorrRes = sum(crossCorrSeq2 .* conj(crossCorrSeq1));
freqOffsetEst = angle(crossCorrRes)/(2*pi)*modRate/interval;

end

