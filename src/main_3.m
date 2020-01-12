clear;
close all;

%% Parameters

% frame parameters
syncLen = 64;
dataLen = 256;
frameLen = 2*syncLen + dataLen;

% physical layer parameters
modRate = 16e6;
clkFreq = 16e6;
sps = clkFreq/modRate;

% channel parameters
phaseOffset = 0;
% freqOffset = 2e3;

norFreq = (-0.6:0.03:0.6)';
EbNo = -4;
repeatTimes = 5000;

GmskMod = comm.GMSKModulator('BitInput', true, 'SamplesPerSymbol', sps, ...
    'PulseLength', 1);

% DFT Algorithm
dftNorFreqOffsetEstTemp = zeros(length(norFreq), repeatTimes);
dftNorFreqOffsetEst = zeros(length(norFreq), 1);

% Kay Algorithm
KayNorFreqOffsetEstTemp = zeros(length(norFreq), repeatTimes);
KayNorFreqOffsetEst = zeros(length(norFreq), 1);

% Fitz Algorithm
FitzNorFreqOffsetEstTemp = zeros(length(norFreq), repeatTimes);
FitzNorFreqOffsetEst = zeros(length(norFreq), 1);

% M&M Algorithm
MaMNorFreqOffsetEstTemp = zeros(length(norFreq), repeatTimes);
MaMNorFreqOffsetEst = zeros(length(norFreq), 1);

% Hybrid Algorithm
HybridNorFreqOffsetEstTemp = zeros(length(norFreq), repeatTimes);
HybridNorFreqOffsetEst = zeros(length(norFreq), 1);

%% Simulation
for i = 1:length(norFreq)
    
    freqOffset = modRate * norFreq(i);
    fprintf('Frequency Offset = %4.0fKHz ...\n', freqOffset/1000);
    channel = comm.AWGNChannel('EbNo', EbNo, 'BitsPerSymbol', 1);
    
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
        spAddNoiseSig = channel(spGmskModSig); % add noise
        spAddPhaseOffsetSig = spAddNoiseSig .* exp(1j*phaseOffset); % add phase offset
        spRxGmskSig = spAddPhaseOffsetSig .* ...
            exp(1j*2*pi*freqOffset*(0:sps*frameLen-1)'/clkFreq);
        
%         dpAddNoiseSig = channel(dpGmskModSig); % add noise
%         dpAddPhaseOffsetSig = dpAddNoiseSig .* exp(1j*phaseOffset); % add phase offset
%         dpRxGmskSig = dpAddPhaseOffsetSig .* ...
%             exp(1j*2*pi*freqOffset*(0:sps*frameLen-1)'/clkFreq);
        
        %% Receiver
        %         decRxGMSK = downsample(rxSigGMSK, sps, sps/2);
        decRxGmskSig = downsample(spRxGmskSig, sps);
        dephaseRx = decRxGmskSig .* conj(spGmskModSig);
%         dephasePre = dephaseRx(1:syncLen);
%         dephasePost = dephaseRx(syncLen+dataLen+1:frameLen);
        
        dephasePilot = dephaseRx(1:syncLen*2);
%         rxSyncPre = decRxGMSK(1:syncLen);
%         rxSyncPost = decRxGMSK(syncLen+dataLen+1:frameLen);
%         
%         localSyncPre = sigGMSKMod(1:syncLen);
%         localSyncPost = sigGMSKMod(syncLen+dataLen+1:frameLen);
%         dephasePre = rxSyncPre .* conj(localSyncPre);
%         dephasePost = rxSyncPost .* conj(localSyncPost);
%         %% DFT Algorithm
%         dftNorFreqOffsetEstTemp(i, time) = dftFreqEstimate(dephasePilot, modRate, 1024) ... 
%             / modRate;
        
        %% Kay Algorithm
        KayNorFreqOffsetEstTemp(i, time) = KayFreqEstimate(dephasePilot, modRate) ...
            / modRate;
        
        %% Fitz Algorithm
        FitzNorFreqOffsetEstTemp(i, time) = FitzFreqEstimate(dephasePilot, modRate) ...
            / modRate;

        %% M&M Algorithm
        MaMNorFreqOffsetEstTemp(i, time) = MaMFreqEstimate(dephasePilot, modRate) ...
            / modRate;
        
        %% New Algorithm
        HybridNorFreqOffsetEstTemp(i, time) = HybridFreqEstimate(dephasePilot, modRate) ...
            / modRate;
    end
    
    
%     dftNorFreqOffsetEst(i) = mean(dftNorFreqOffsetEstTemp(i, :));
    KayNorFreqOffsetEst(i) = mean(KayNorFreqOffsetEstTemp(i, :));
    MaMNorFreqOffsetEst(i) = mean(MaMNorFreqOffsetEstTemp(i, :));
    FitzNorFreqOffsetEst(i) = mean(FitzNorFreqOffsetEstTemp(i, :));
    HybridNorFreqOffsetEst(i) = mean(HybridNorFreqOffsetEstTemp(i, :));
    
end

%% Performance Compare

dftNorFreqVar = zeros(length(norFreq), 1);
KayNorFreqVar = zeros(length(norFreq), 1);
MaMNorFreqVar = zeros(length(norFreq), 1);
FitzNorFreqVar = zeros(length(norFreq), 1);
HybridNorFreqVar = zeros(length(norFreq), 1);

for i = 1:length(norFreq)
    
    dftNorFreqErr = dftNorFreqOffsetEstTemp(i, :) - norFreq(i);
    KayNorFreqErr = KayNorFreqOffsetEstTemp(i, :) - norFreq(i);
    MaMNorFreqErr = MaMNorFreqOffsetEstTemp(i, :) - norFreq(i);
    FitzNorFreqErr = FitzNorFreqOffsetEstTemp(i, :) - norFreq(i);
    HybridNorFreqErr = HybridNorFreqOffsetEstTemp(i, :) - norFreq(i);
    
    dftNorFreqVar(i) = sum(dftNorFreqErr .^ 2)/repeatTimes;
    KayNorFreqVar(i) = sum(KayNorFreqErr .^ 2)/repeatTimes;
    MaMNorFreqVar(i) = sum(MaMNorFreqErr .^ 2)/repeatTimes;
    FitzNorFreqVar(i) = sum(FitzNorFreqErr .^ 2)/repeatTimes;
    HybridNorFreqVar(i) = sum(HybridNorFreqErr .^ 2)/repeatTimes;
end

figure;
plot(norFreq, KayNorFreqOffsetEst, '-s'); hold on
plot(norFreq, MaMNorFreqOffsetEst, '-x'); hold on
plot(norFreq, FitzNorFreqOffsetEst, '-*'); hold on
plot(norFreq, HybridNorFreqOffsetEst, '-o'); hold on
legend('Kay', 'M&M', 'Fitz', 'Hybrid');
title(['Frequency Offset Scale(EbNo = ', num2str(EbNo), 'dB)']);

figure;
semilogy(norFreq, KayNorFreqVar, '-s'); hold on
semilogy(norFreq, MaMNorFreqVar, '-x'); hold on
semilogy(norFreq, FitzNorFreqVar, '-*'); hold on
semilogy(norFreq, HybridNorFreqVar, '-o'); hold on
legend('Kay', 'M&M', 'Fitz', 'Hybrid');
title(['VAR of Normalization Frequency Offset(EbNo = ', num2str(EbNo), 'dB)']);
% savefig(['Frequency Offset Scale(EbNo=', num2str(EbNo), 'dB).fig']);