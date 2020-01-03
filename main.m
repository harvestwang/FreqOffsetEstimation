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
freqOffset = 1e3;
EbNo = (-15:15)';

alpha = 1;
repeatTimes = 3000;

firFreqOffsetEst = zeros(length(EbNo), 1);
secFreqOffsetEst = zeros(length(EbNo), 1);

% newFreqOffsetEst = zeros(length(EbNo), 1);

gmskMod = comm.GMSKModulator('BitInput', true, 'SamplesPerSymbol', sps, ...
                                'PulseLength', 1);

%% Simulation
for i = 1:length(EbNo)

    fprintf("\nEbNo = %d ...", EbNo(i));
    
    preFreqOffsetEst = zeros(repeatTimes, 1);
    secFreqOffsetEstTemp = zeros(repeatTimes, 1);
    
    for time = 1 : repeatTimes
    %% Initialization
        [syncPreSrc, dataSrc, syncPostSrc, syncPreCode, ...
            dataCode, syncPostCode] = sourceGen(syncLen, dataLen);

    %% Transmitter
        frame = [syncPreCode; dataCode; syncPostCode];
        sigGMSKMod = gmskMod(frame);
%         txSigGMSK = downsample(sigGMSKMod, sps, sps/2);
%         scatterplot(sigGMSKMod);
    %% Channel
        channel = comm.AWGNChannel('EbNo', EbNo(i), 'BitsPerSymbol', 1);
        addNoiseSig = channel(sigGMSKMod); % add noise
        addPhaseOffsetSig = addNoiseSig .* exp(1j*phaseOffset); % add phase offset
        rxSigGMSK = addPhaseOffsetSig .* ...
            exp(1j*2*pi*freqOffset*(0:sps*frameLen-1)'/clkFreq);
    
    %% Receiver
%         decRxGMSK = downsample(rxSigGMSK, sps, sps/2);
        decRxGMSK = downsample(rxSigGMSK, sps);
%         scatterplot(decRxGMSK);
        decLocalGMSK = sigGMSKMod;
        
        rxSignal = decRxGMSK .* conj(sigGMSKMod);
        
        % Fisrt frequency offset estimation
        
        % L&W Alogithm
        k = (1 : syncLen-1)';
        w1 = 6 .* k .* (syncLen-k)/(syncLen*(syncLen-1)); % window
        
        selfCorrPreSeq1 = rxSignal(1 : syncLen-alpha);
        selfCorrPreSeq2 = rxSignal(alpha+1 : syncLen);
        selfCorrPreRes = sum(w1 .* selfCorrPreSeq2 .* conj(selfCorrPreSeq1));
        preFreqOffsetEst(time) = angle(selfCorrPreRes)/(2*pi)*modRate/alpha;
        
%         selfCorrPostSeq1 = rxSignal(syncLen+dataLen+1 : frameLen-alpha);
%         selfCorrPostSeq2 = rxSignal(alpha+syncLen+dataLen+1 : frameLen);
%         selfCorrPostRes = sum(w1 .* selfCorrPostSeq2 .* conj(selfCorrPostSeq1));
%         postFreqOffsetEst = angle(selfCorrPostRes)/(2*pi)*modRate/alpha;
        
%         firFreqOffsetEst = preFreqOffsetEst;
       
    end
    
%      firFreqOffsetEst = (preFreqOffsetEst + postFreqOffsetEst)/2;
    firFreqOffsetEst(i) = mean(preFreqOffsetEst);
    
    % Fisrt frequency offset correction
    decRxSigCorrect1 = rxSignal .* ...
        exp(-1j*2*pi*(0:frameLen-1)'*firFreqOffsetEst(i)/clkFreq);
%         scatterplot(decRxSigCorrect1);
    
    for time = 1 : repeatTimes
        % Second frequency offset estimation
        
        crossCorrSeq1 = decRxSigCorrect1(1 : syncLen);
        crossCorrSeq2 = decRxSigCorrect1(syncLen+dataLen+1 : frameLen);
        
        crossCorrRes = sum(crossCorrSeq2 .* conj(crossCorrSeq1));
        secFreqOffsetEstTemp(time) = angle(crossCorrRes)/(2*pi)*modRate/(syncLen + dataLen);
    end
    
    secFreqOffsetEst(i) = mean(secFreqOffsetEstTemp);
    
    % second frequency offset correction
    decRxSigCorrect2 = decRxSigCorrect1 .* ... 
        exp(-1j*2*pi*(0:frameLen-1)'*secFreqOffsetEst(i)/clkFreq);
%         scatterplot(decRxSigCorrect2);
end

newFreqOffsetEst = firFreqOffsetEst + secFreqOffsetEst;
%     norFreqOffset = mean((estFreqOffsetTemp - freqOffset).^2);
%     newFreqOffsetEst(i) = sqrt(norFreqOffset)/freqOffset;

%% Plot results
fprintf("\n");

figure;
plot(EbNo, newFreqOffsetEst, '-x'); hold on
plot(EbNo, firFreqOffsetEst, '-s'); hold on
plot(EbNo, secFreqOffsetEst, '-*'); hold on
% plot(EbNo, estEbNoPaper_dB, '-o'); hold on
% plot(EbNo, estEbNoUnbias_dB, '-^'); hold on
% plot(EbNo, preEstEbNoUnbias_dB, '-<');
% plot(EbNo, estEbNoSelfCorr_dB, '->');
% legend('Ideal', 'GMSK Estimation', 'TTNT Estimation', 'Paper Estimation', ...
%         'Unbias Estimation', 'Pre Unbias Estimation');
% title("SNR estimation");
% 
% estDeltaGMSK = abs(estEbNoGMSK_dB - EbNo);
% estDeltaTTNT = abs(estEbNoTTNT_dB - EbNo);
% estDeltaPaper = abs(estEbNoPaper_dB - EbNo);
% estDeltaUnbias = abs(estEbNoUnbias_dB - EbNo);
% preEstDeltaUnbias = abs(preEstEbNoUnbias_dB - EbNo);
% 
% figure;
% plot(EbNo, estDeltaGMSK, '-s'); hold on
% plot(EbNo, estDeltaTTNT, '-*'); hold on
% plot(EbNo, estDeltaPaper, '-o'); hold on
% plot(EbNo, estDeltaUnbias, '-^'); hold on
% plot(EbNo, preEstDeltaUnbias, '-<');
% legend('GMSK Error', 'TTNT Error', 'Paper Error', 'Unbiased Error', ...
%         'Pre Unbiased Error');
% title("Estimated SNR Error");