function [fftRes,dftFreqEst] = dftFreqEstimate(inSignal, Fs, fftN)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    fftRes = fftshift(fft(inSignal, fftN));
    absFftRes = abs(fftRes);
    fftIndex = find(absFftRes == max(fftRes));
    fprintf("fftIndex = %d\n", fftIndex);
    dftFreqEst = (fftIndex-fftN/2)/fftN*Fs;
end