function freqOffsetEst = FitzFreqEstimate(pilot, modRate)

    pilotLen = length(pilot);
    N = pilotLen/2;
    phaseDelta = zeros(N, 1);
    
    for i = 1:N
        phaseDelta(i) = angle(sum(pilot(i+1:pilotLen) .* conj(pilot(1:pilotLen-i))) ...
            / (pilotLen-i));
    end
    
    freqOffsetEst = sum(phaseDelta)*modRate/(pi*N*(N+1));
end

