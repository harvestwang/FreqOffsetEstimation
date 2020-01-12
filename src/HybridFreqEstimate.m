function freqOffsetEst = HybridFreqEstimate(pilot, modRate)
    
    pilotLen = length(pilot);
    fftN = 2^(ceil(log2(pilotLen))+1);

    firFreqOffsetEst = dftFreqEstimate(pilot, modRate, fftN);
    
%     fftRes = fftshift(fft(pilot, fftN));
%     absFftRes = abs(fftRes);
%     fftIndex = find(absFftRes == max(absFftRes));
%     firFreqOffsetEst = (fftIndex-fftN/2)/fftN*Fs;
    
    correctedPilot = pilot .* ...
        exp(-1j*2*pi*(0:pilotLen-1)'*firFreqOffsetEst/modRate);
    
%     secFreqOffsetEst = FitzFreqEstimate(correctedPilot, modRate);
    N = pilotLen/2;
    
    phaseDelta = zeros(N, 1);
    for i = 1:N
        phaseDelta(i) = angle(sum(correctedPilot(i+1:pilotLen) .* ...
            conj(correctedPilot(1:pilotLen-i))));
    end
    secFreqOffsetEst = sum(phaseDelta)*modRate/(pi*N*(N+1));
    
    freqOffsetEst = firFreqOffsetEst + secFreqOffsetEst;
end

