function freqOffsetEst = KayFreqEstimate(pilot, modRate)
    
    pilotLen = length(pilot);
    k = (1 : pilotLen-1)';
    w2 = (3/2)*(pilotLen/(pilotLen^2-1))*(1-(2*k/pilotLen-1).^2);

    phaseDelta = angle(pilot(2:pilotLen) .* conj(pilot(1:pilotLen-1)));
    
    freqOffsetEst = sum(w2 .* phaseDelta)*modRate/(2*pi);
end

