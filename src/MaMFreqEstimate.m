function freqOffsetEst = MaMFreqEstimate(pilot, modRate)

    pilotLen = length(pilot);
    N = pilotLen/2;
    m = (1:N)';
    w = 3*((pilotLen-m).*(pilotLen-m+1)-N*(pilotLen-N))./N/(4*N^2-6*N*pilotLen+3*pilotLen^2-1);


    phaseDelta = zeros(N, 1);
    phaseDelta(1) = angle(sum(pilot(2:pilotLen) .* conj(pilot(1:pilotLen-1))));
    
    for i = 2:N
    	
    	R1 = sum(pilot(i+1:pilotLen) .* conj(pilot(1:pilotLen-i)))/(pilotLen-i);
    	R2 = sum(pilot(i:pilotLen) .* conj(pilot(1:pilotLen-i+1)))/(pilotLen-i-1);
    	
        if (angle(R1) - angle(R2)) < -pi
            phaseDelta(i) = angle(R1)-angle(R2)+2*pi;
        elseif (angle(R1) - angle(R2)) > pi
            phaseDelta(i) = angle(R1)-angle(R2)-2*pi;
        else
            phaseDelta(i) = angle(R1)-angle(R2);
        end
    end

    freqOffsetEst = sum(w .* phaseDelta)*modRate/(2*pi);

end

