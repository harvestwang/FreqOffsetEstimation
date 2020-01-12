clear;
close all;

%% Parameters

% frame parameters
dataLen = 256;
syncLen = (10:1:40)';
spLen = 2*syncLen;

% physical layer parameters
modRate = 16e6;
clkFreq = 16e6;
sps = clkFreq/modRate;

% channel parameters
phaseOffset = 0;
% freqOffset = 2e3;
norFreq = 0.005;
freqOffset = modRate * norFreq;

EbNo = 10;
repeatTimes = 10000;

% GmskMod = comm.GMSKModulator('BitInput', true, 'SamplesPerSymbol', sps, ...
%     'PulseLength', 1);


% DFT Algorithm
dftNorFreqOffsetEstTemp = zeros(length(syncLen), repeatTimes);
dftNorFreqOffsetEst = zeros(length(syncLen), 1);

% Kay Algorithm
KayNorFreqOffsetEstTemp = zeros(length(syncLen), repeatTimes);
KayNorFreqOffsetEst = zeros(length(syncLen), 1);

% Fitz Algorithm
FitzNorFreqOffsetEstTemp = zeros(length(syncLen), repeatTimes);
FitzNorFreqOffsetEst = zeros(length(syncLen), 1);

% M&M Algorithm
MaMNorFreqOffsetEstTemp = zeros(length(syncLen), repeatTimes);
MaMNorFreqOffsetEst = zeros(length(syncLen), 1);

% Hybrid Algorithm
HybridNorFreqOffsetEstTemp = zeros(length(syncLen), repeatTimes);
HybridNorFreqOffsetEst = zeros(length(syncLen), 1);

%% Simulation
for i = 1:length(spLen)
    
    frameLen = spLen(i) + dataLen;
    
    fprintf('Pilot Length = %d.\n', spLen(i));
    
    GmskMod = comm.GMSKModulator('BitInput', true, 'SamplesPerSymbol', sps, ...
        'PulseLength', 1);
    
    for time = 1 : repeatTimes
        %% Initialization
        [syncPreSrc, dataSrc, syncPostSrc, syncPreCode, ...
            dataCode, syncPostCode] = sourceGen(syncLen(i), dataLen);
        
        %% Transmitter
        spFrame = [syncPreCode; syncPostCode; dataCode]; % single pilot L=48
%         dpFrame = [syncPreCode; dataCode; syncPostCode]; % double pilot Lpre=Lpost=24
        
        spGmskModSig = GmskMod(spFrame);
%         dpGmskModSig = GmskMod(dpFrame);
        
        %% Channel
        channel = comm.AWGNChannel('EbNo', EbNo, 'BitsPerSymbol', 1);
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
        
        dephasePilot = dephaseRx(1:spLen(i));
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

dftNorFreqVar = zeros(length(syncLen), 1);
KayNorFreqVar = zeros(length(syncLen), 1);
MaMNorFreqVar = zeros(length(syncLen), 1);
FitzNorFreqVar = zeros(length(syncLen), 1);
HybridNorFreqVar = zeros(length(syncLen), 1);

for i = 1:length(syncLen)
    
    dftNorFreqErr = dftNorFreqOffsetEstTemp(i, :) - norFreq;
    KayNorFreqErr = KayNorFreqOffsetEstTemp(i, :) - norFreq;
    MaMNorFreqErr = MaMNorFreqOffsetEstTemp(i, :) - norFreq;
    FitzNorFreqErr = FitzNorFreqOffsetEstTemp(i, :) - norFreq;
    HybridNorFreqErr = HybridNorFreqOffsetEstTemp(i, :) - norFreq;
    
    dftNorFreqVar(i) = sum(dftNorFreqErr .^ 2)/repeatTimes;
    KayNorFreqVar(i) = sum(KayNorFreqErr .^ 2)/repeatTimes;
    MaMNorFreqVar(i) = sum(MaMNorFreqErr .^ 2)/repeatTimes;
    FitzNorFreqVar(i) = sum(FitzNorFreqErr .^ 2)/repeatTimes;
    HybridNorFreqVar(i) = sum(HybridNorFreqErr .^ 2)/repeatTimes;
end

% figure;
% plot(norFreq, KayNorFreqOffsetEst, '-s'); hold on
% plot(norFreq, MaMNorFreqOffsetEst, '-x'); hold on
% plot(norFreq, FitzNorFreqOffsetEst, '-*'); hold on
% plot(norFreq, HybridNorFreqOffsetEst, '-o'); hold on
% legend('Kay', 'M&M', 'Fitz', 'Hybrid');
% title(['Frequency Offset Scale(PilotLen = ', num2str(2*syncLen), 'dB)']);

decEbNo = 10^(EbNo/10);
CRB = 6/(2*pi)^2./spLen./(spLen.^2-1)/decEbNo;

figure;
semilogy(spLen, KayNorFreqVar, '-s'); hold on
semilogy(spLen, MaMNorFreqVar, '-x'); hold on
semilogy(spLen, FitzNorFreqVar, '-<'); hold on
semilogy(spLen, HybridNorFreqVar, '-o'); hold on
semilogy(spLen, CRB, '-*'); hold on
legend('Kay', 'M&M', 'Fitz', 'Hybrid', 'CRB');
title('VAR of Normalization Frequency Offset');
% savefig(['Frequency Offset Scale(EbNo=', num2str(EbNo), 'dB).fig']);