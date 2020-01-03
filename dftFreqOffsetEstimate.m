function [fftRes,dftFreqEst] = dftFreqOffsetEstimate(inSignal, Fs, fftN)

    fftRes = fftshift(fft(inSignal, fftN));
    absFftRes = abs(fftRes);
    fftIndex = find(absFftRes == max(absFftRes));
%     fprintf('fftIndex = %d\n', fftIndex);
    dftFreqEst = (fftIndex-fftN/2)/fftN*Fs;
end