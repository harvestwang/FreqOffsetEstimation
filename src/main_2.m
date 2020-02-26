clear;
close all;

%% Parameters

% frame parameters
syncLen = 32;
dataLen = 256;
spLen = 2*syncLen;
frameLen = 2*syncLen + dataLen;

% physical layer parameters
modRate = 16e6;
clkFreq = 16e6;
sps = clkFreq/modRate;

% channel parameters
norFreqOffset = 0.35;
phaseOffset = 0;
freqOffset = norFreqOffset*modRate;
SNR = (-10:10)';

repeatTimes = 10000;

GmskMod = comm.GMSKModulator('BitInput', true, 'SamplesPerSymbol', sps, ...
    'PulseLength', 1);

% Hybrid Algorithm
HybridFreqOffsetEstTemp = zeros(length(SNR), repeatTimes);
HybridFreqOffsetEst = zeros(length(SNR), 1);

% New Algorithm
newFreqOffsetEstTemp = zeros(length(SNR), repeatTimes);
newFreqOffsetEst = zeros(length(SNR), 1);

% DFT Algorithm
dftFreqOffsetEstTemp = zeros(length(SNR), repeatTimes);
dftFreqOffsetEst = zeros(length(SNR), 1);

% Kay Algorithm
KayFreqOffsetEstTemp = zeros(length(SNR), repeatTimes);
KayFreqOffsetEst = zeros(length(SNR), 1);

% M&M Algorithm
MaMFreqOffsetEstTemp = zeros(length(SNR), repeatTimes);
MaMFreqOffsetEst = zeros(length(SNR), 1);

% Fitz Algorithm
FitzFreqOffsetEstTemp = zeros(length(SNR), repeatTimes);
FitzFreqOffsetEst = zeros(length(SNR), 1);

%% Simulation
for i = 1:length(SNR)
    
    fprintf('EbNo = %2ddB ...\n', SNR(i));
%     channel = comm.AWGNChannel('EbNo', EbNo(i), 'BitsPerSymbol', 1);
    
    for time = 1 : repeatTimes
        %% Initialization
        [syncPreSrc, dataSrc, syncPostSrc, syncPreCode, ...
            dataCode, syncPostCode] = sourceGen(syncLen, dataLen);
        
        %% Transmitter
        spFrame = [syncPreCode; syncPostCode; dataCode]; % single pilot L=48
%         dpFrame = [syncPreCode; dataCode; syncPostCode]; % double pilot Lpre=Lpost=24
        
        spGmskModSig = GmskMod(spFrame);
%         dpGmskModSig = GmskMod(dpFrame);
        
        %% Channel
%         spAddNoiseSig = channel(spGmskModSig); % add noise
        spAddNoiseSig = awgn(spGmskModSig, SNR(i));
        spAddPhaseOffsetSig = spAddNoiseSig .* exp(1j*phaseOffset); % add phase offset
        spRxGmskSig = spAddPhaseOffsetSig .* ...
            exp(1j*2*pi*freqOffset*(0:sps*frameLen-1)'/clkFreq);
        
%         dpAddNoiseSig = channel(dpGmskModSig); % add noise
%         dpAddPhaseOffsetSig = dpAddNoiseSig .* exp(1j*phaseOffset); % add phase offset
%         dpRxGmskSig = dpAddPhaseOffsetSig .* ...
%             exp(1j*2*pi*freqOffset*(0:sps*frameLen-1)'/clkFreq);
        
        %% Receiver
        % Single pilot frame
        spRxGmskSig = downsample(spRxGmskSig, sps);
        spDephaseRx = spRxGmskSig .* conj(spGmskModSig);
        spDephasePilot = spDephaseRx(1:spLen);
        
        % Double pilot frame 
%         dpRxGmskSig = downsample(dpRxGmskSig, sps);
%         dpDephaseRx = dpRxGmskSig .* conj(dpGmskModSig);
        
%         dpDephasePre = dpDephaseRx(1:syncLen);
%         dpDephasePost = dpDephaseRx(syncLen+dataLen+1:frameLen);

        %% DFT Algorithm
        dftFreqOffsetEstTemp(i, time) = dftFreqEstimate(spDephasePilot, modRate, 1024);
        
        %% Kay Algorithm
%         KayFreqOffsetEstTemp(i, time) = KayFreqEstimate(spDephasePilot, modRate);
        
        %% M&M Algorithm
        MaMFreqOffsetEstTemp(i, time) = MaMFreqEstimate(spDephasePilot, modRate);
        
        %% Fitz Algorithm
        FitzFreqOffsetEstTemp(i, time) = FitzFreqEstimate(spDephasePilot, modRate);
        
        %% Hybrid Algorithm
        HybridFreqOffsetEstTemp(i, time) = HybridFreqEstimate(spDephasePilot, modRate);

    end
    
    dftFreqOffsetEst(i) = mean(dftFreqOffsetEstTemp(i, :));
%     KayFreqOffsetEst(i) = mean(KayFreqOffsetEstTemp(i, :));
    MaMFreqOffsetEst(i) = mean(MaMFreqOffsetEstTemp(i, :));
    FitzFreqOffsetEst(i) = mean(FitzFreqOffsetEstTemp(i, :));
    HybridFreqOffsetEst(i) = mean(HybridFreqOffsetEstTemp(i, :));
    
end

%% Performance Compare

% NRMSE
dftNorRMSE = zeros(length(SNR), 1);
% KayNorRMSE = zeros(length(EbNo), 1);
MaMNorRMSE = zeros(length(SNR), 1);
FitzNorRMSE = zeros(length(SNR), 1);
HybridNorRMSE = zeros(length(SNR), 1);

for i = 1:length(SNR)
    dftNorRMSE(i) = sqrt(mean((dftFreqOffsetEstTemp(i,:)-freqOffset).^2))/freqOffset;
%     KayNorRMSE(i) = sqrt(mean((KayFreqOffsetEstTemp(i,:)-freqOffset).^2))/freqOffset;
    MaMNorRMSE(i) = sqrt(mean((MaMFreqOffsetEstTemp(i,:)-freqOffset).^2))/freqOffset;
    FitzNorRMSE(i) = sqrt(mean((FitzFreqOffsetEstTemp(i,:)-freqOffset).^2))/freqOffset;
    HybridNorRMSE(i) = sqrt(mean((HybridFreqOffsetEstTemp(i,:)-freqOffset).^2))/freqOffset;
end

% NVAR
dftNorFreqVar = zeros(length(SNR), 1);
% KayNorFreqVar = zeros(length(EbNo), 1);
MaMNorFreqVar = zeros(length(SNR), 1);
FitzNorFreqVar = zeros(length(SNR), 1);
HybridNorFreqVar = zeros(length(SNR), 1);
for i = 1:length(SNR)
    
    dftNorFreq = (dftFreqOffsetEstTemp(i,:)-freqOffset)/modRate;
%     KayNorFreq = (KayFreqOffsetEstTemp(i,:)-freqOffset)/modRate;
    MaMNorFreq = (MaMFreqOffsetEstTemp(i,:)-freqOffset)/modRate;
    FitzNorFreq = (FitzFreqOffsetEstTemp(i,:)-freqOffset)/modRate;
    HybridNorFreq = (HybridFreqOffsetEstTemp(i,:)-freqOffset)/modRate;
    
    dftNorFreqVar(i) = sum(dftNorFreq .^ 2)/repeatTimes;
%     KayNorFreqVar(i) = sum(KayNorFreq .^ 2)/repeatTimes;
    MaMNorFreqVar(i) = sum(MaMNorFreq .^ 2)/repeatTimes;
    FitzNorFreqVar(i) = sum(FitzNorFreq .^ 2)/repeatTimes;
    HybridNorFreqVar(i) = sum(HybridNorFreq .^ 2)/repeatTimes;
end

decSNR = 10.^(SNR/10);
CRB = 6/(2*pi)^2/spLen/(spLen^2-1)./decSNR;
% figure;semilogy(EbNo, CRB);
% save(['newNorRMSE_', num2str(freqOffset/1000), 'kHz.mat'], 'newNorRMSE');
% save(['newNorFreqVar_', num2str(freqOffset/1000), 'kHz.mat'], 'newNorFreqVar');
% 
% save(['dftNorRMSE_', num2str(freqOffset/1000), 'kHz.mat'], 'dftNorRMSE');
% save(['dftNorFreqVar_', num2str(freqOffset/1000), 'kHz.mat'], 'dftNorFreqVar');
% 
% save(['KayNorRMSE_', num2str(freqOffset/1000), 'kHz.mat'], 'KayNorRMSE');
% save(['KayNorFreqVar_', num2str(freqOffset/1000), 'kHz.mat'], 'KayNorFreqVar');
% 
% save(['FitzNorRMSE_', num2str(freqOffset/1000), 'kHz.mat'], 'FitzNorRMSE');
% save(['FitzNorFreqVar_', num2str(freqOffset/1000), 'kHz.mat'], 'FitzNorFreqVar');

% figure;
% plot(EbNo, newFreqOffsetEst, '-o'); hold on
% plot(EbNo, dftFreqOffsetEst, '-x'); hold on
% plot(EbNo, KayFreqOffsetEst, '-s'); hold on
% plot(EbNo, FitzFreqOffsetEst, '-*');
% legend('New', 'DFT', 'Kay', 'Fitz');
% title(['Frequency Offset Esitamation(offset = ', num2str(freqOffset/1000), 'kHz)']);
% savefig(['Frequency(offset = ', num2str(freqOffset/1000), 'kHz).fig']);

% figure;
% plot(EbNo, KayNorRMSE, '-s'); hold on
% plot(EbNo, MaMNorRMSE, '-x'); hold on
% plot(EbNo, FitzNorRMSE, '-<'); hold on
% plot(EbNo, HybridNorRMSE, '-o'); hold on
% legend('Kay', 'M&M', 'Fitz', 'Hybrid');
% title(['NRMSE(offset = ', num2str(freqOffset/1000), 'kHz)']);
% savefig(['NRMSE(offset = ', num2str(freqOffset/1000), 'kHz).fig']);

figure;
semilogy(SNR, dftNorFreqVar, '-s'); hold on
semilogy(SNR, MaMNorFreqVar, '-x'); hold on
semilogy(SNR, FitzNorFreqVar, '-<'); hold on
semilogy(SNR, HybridNorFreqVar, '-o'); hold on
semilogy(SNR, CRB, '-*'); hold on
legend('FFT(N_{fft}=1024)', 'M&M', 'Fitz', 'proposed', 'CRB');
xlabel('SNR/dB');
ylabel('var(f_dT_s)');
% title(['f_dT_s = ', num2str(freqOffset/modRate)]);
% savefig(['..\', 'NVAR(offset = ', num2str(freqOffset/1000), 'kHz).fig']);