clear;
close all;

%% Parameters

% frame parameters
syncLen = 24;
dataLen = 256;
frameLen = 2*syncLen + dataLen;

% physical layer parameters
modRate = 16e6;
clkFreq = 16e6;
sps = clkFreq/modRate;

% channel parameters
phaseOffset = 0;
% freqOffset = 2e3;

norFreq = (-0.02:0.002:0.02)';
EbNo = 10;
repeatTimes = 5000;

GmskMod = comm.GMSKModulator('BitInput', true, 'SamplesPerSymbol', sps, ...
    'PulseLength', 1);

% new Algorithm
% firFreqOffsetEstTemp = zeros(length(EbNo), repeatTimes);
% secFreqOffsetEstTemp = zeros(length(EbNo), repeatTimes);
newNorFreqOffsetEstTemp = zeros(length(EbNo), repeatTimes);
newNorFreqOffsetEst = zeros(length(EbNo), 1);

% secFreqOffsetEstTemp1 = zeros(length(EbNo), repeatTimes);
% newFreqOffsetEstTemp1 = zeros(length(EbNo), repeatTimes);
% newFreqOffsetEst1 = zeros(length(EbNo), 1);

% DFT Algorithm
dftNorFreqOffsetEstTemp = zeros(length(EbNo), repeatTimes);
dftNorFreqOffsetEst = zeros(length(EbNo), 1);

% Kay Algorithm
KayNorFreqOffsetEstTemp = zeros(length(EbNo), repeatTimes);
KayNorFreqOffsetEst = zeros(length(EbNo), 1);

% Fitz Algorithm
FitzNorFreqOffsetEstTemp = zeros(length(EbNo), repeatTimes);
FitzNorFreqOffsetEst = zeros(length(EbNo), 1);
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
        %% DFT Algorithm
        dftNorFreqOffsetEstTemp(i, time) = dftFreqEstimate(dephasePilot, modRate, 1024) ... 
            / modRate;
        
        %% Kay Algorithm
        KayNorFreqOffsetEstTemp(i, time) = KayFreqEstimate(dephasePilot, modRate) ...
            / modRate;
        %% Fitz Algorithm
        FitzNorFreqOffsetEstTemp(i, time) = FitzFreqEstimate(dephasePilot, modRate) ...
            / modRate;
        %% New Algorithm
        newNorFreqOffsetEstTemp(i, time) = newFreqEstimate(dephasePilot, modRate) ...
            / modRate;
        % Fisrt frequency offset estimation and correction
        
%         k = (1 : syncLen-1)';
%         w1 = 6 .* k .* (syncLen-k)/(syncLen*(syncLen-1)); % window
%         
%         selfCorrPreSeq1 = dephasePre(1 : syncLen-alpha);
%         selfCorrPreSeq2 = dephasePre(alpha+1 : syncLen);
%         selfCorrPreRes = mean(w1 .* selfCorrPreSeq2 .* conj(selfCorrPreSeq1));
%         firFreqOffsetEstTemp(i, time) = angle(selfCorrPreRes)/(2*pi)*modRate/alpha;
        
%         dephaseRxCorrected = dephaseRx .* ...
%             exp(-1j*2*pi*(0:frameLen-1)'*firFreqOffsetEstTemp(i, time)/clkFreq);
%         rxSyncPreCorrect = dephasePre .* ...
%             exp(-1j*2*pi*(0:syncLen-1)'*firFreqOffsetEstTemp(i, time)/clkFreq);
%         rxSyncPostCorrect = dephasePost .* ...
%             exp(-1j*2*pi*(syncLen+dataLen:frameLen-1)'*firFreqOffsetEstTemp(i, time)/clkFreq);
        
        % Second frequency offset estimation and correction
%         crossCorrSeq1 = dephaseRxCorrected(1 : syncLen);
%         crossCorrSeq2 = dephaseRxCorrected(syncLen+dataLen+1 : frameLen);
%         interval = syncLen + dataLen;
        
%         secFreqOffsetEstTemp(i, time) = crossCorrFreqEstimate(crossCorrSeq1, ...
%             crossCorrSeq2, modRate, interval);
        
%         dephaseRxCorrected1 = dephaseRx .* ...
%             exp(-1j*2*pi*(0:frameLen-1)'*dftFreqOffsetEstTemp(i, time)/clkFreq);
        
%         crossCorrSeq1 = dephaseRxCorrected1(1 : syncLen);
%         crossCorrSeq2 = dephaseRxCorrected1(syncLen+dataLen+1 : frameLen);
%         interval = syncLen + dataLen;
        
%         secFreqOffsetEstTemp1(i, time) = crossCorrFreqEstimate(crossCorrSeq1, ...
%             crossCorrSeq2, modRate, interval);
        
%         newFreqOffsetEstTemp(i, time) = firFreqOffsetEstTemp(i, time) + secFreqOffsetEstTemp(i, time);
%         newFreqOffsetEstTemp1(i, time) = dftFreqOffsetEstTemp(i, time) + secFreqOffsetEstTemp1(i, time);

        % second frequency offset correction
        %         decRxSigCorrect2 = decRxSigCorrect1 .* ...
        %             exp(-1j*2*pi*(0:syncLen-1)'*secFreqOffsetEstTemp(i, time)/clkFreq);
    end
    
    newNorFreqOffsetEst(i) = mean(newNorFreqOffsetEstTemp(i, :));
    dftNorFreqOffsetEst(i) = mean(dftNorFreqOffsetEstTemp(i, :));
    KayNorFreqOffsetEst(i) = mean(KayNorFreqOffsetEstTemp(i, :));
    FitzNorFreqOffsetEst(i) = mean(FitzNorFreqOffsetEstTemp(i, :));
    
end

%% Performance Compare

newNorFreqVar = zeros(length(norFreq), 1);
dftNorFreqVar = zeros(length(norFreq), 1);
KayNorFreqVar = zeros(length(norFreq), 1);
FitzNorFreqVar = zeros(length(norFreq), 1);
for i = 1:length(norFreq)
    
    newNorFreqErr = newNorFreqOffsetEstTemp(i, :) - norFreq(i);
    dftNorFreqErr = dftNorFreqOffsetEstTemp(i, :) - norFreq(i);
    KayNorFreqErr = KayNorFreqOffsetEstTemp(i, :) - norFreq(i);
    FitzNorFreqErr = FitzNorFreqOffsetEstTemp(i, :) - norFreq(i);
    
    newNorFreqVar(i) = sum(newNorFreqErr .^ 2)/repeatTimes;
    dftNorFreqVar(i) = sum(dftNorFreqErr .^ 2)/repeatTimes;
    KayNorFreqVar(i) = sum(KayNorFreqErr .^ 2)/repeatTimes;
    FitzNorFreqVar(i) = sum(FitzNorFreqErr .^ 2)/repeatTimes;
end

figure;
plot(norFreq, newNorFreqOffsetEst, '-o'); hold on
% plot(norFreq, dftNorFreqOffsetEst, '-x'); hold on
plot(norFreq, KayNorFreqOffsetEst, '-s'); hold on
plot(norFreq, FitzNorFreqOffsetEst, '-*');
legend('New', 'Kay', 'Fitz');
% legend('New', 'DFT', 'Kay', 'Fitz');
title(['Frequency Offset Scale(EbNo = ', num2str(EbNo), 'dB)']);

figure;
semilogy(norFreq, newNorFreqVar, '-o'); hold on
% semilogy(norFreq, dftNorFreqVar, '-x'); hold on
semilogy(norFreq, KayNorFreqVar, '-s'); hold on
semilogy(norFreq, FitzNorFreqVar, '-*');
legend('New', 'Kay', 'Fitz');
% legend('New', 'DFT', 'Kay', 'Fitz');
title(['VAR of Normalization Frequency Offset(EbNo = ', num2str(EbNo), 'dB)']);
savefig(['Frequency Offset Scale(EbNo=', num2str(EbNo), 'dB).fig']);