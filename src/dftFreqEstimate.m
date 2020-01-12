function freqOffsetEst = dftFreqEstimate(pilot, Fs, fftN)

    fftRes = fftshift(fft(pilot, fftN));
    absFftRes = abs(fftRes);
    fftIndex = find(absFftRes == max(absFftRes));
    freqOffsetEst = (fftIndex-fftN/2)*Fs/fftN;
    
end