function [fftRes,dftFreqEst] = dftFreqEstimate(inSignal, Fs, fftN)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    fftRes = fftshift(fft(inSignal, fftN));
    absFftRes = abs(fftRes);
    fftIndex = find(absFftRes == max(fftRes));
    fprintf("fftIndex = %d\n", fftIndex);
    dftFreqEst = (fftIndex-fftN/2)/fftN*Fs;
end