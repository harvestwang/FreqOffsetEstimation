function [freqOffsetEst] = selfCorrFreqEstimate(inSignal, modRate, interval)

    % L&W Algorithm
    len = length(inSignal);
    k = (1 : len-1)';
    w1 = 6 .* k .* (len-k)/(len*(len-1)); % window

    selfCorrPreSeq1 = inSignal(1 : len-interval);
    selfCorrPreSeq2 = inSignal(interval+1 : len);
    selfCorrPreRes = sum(w1 .* selfCorrPreSeq2 .* conj(selfCorrPreSeq1));
    freqOffsetEst = angle(selfCorrPreRes)/(2*pi)*modRate/interval;
end
