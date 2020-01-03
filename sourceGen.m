  function [syncPreSrc, dataSrc, syncPostSrc, ...
      syncPreCode, dataCode, syncPostCode] = sourceGen(syncLen, dataLen)

    frameLen = syncLen*2 + dataLen;
  
%% 源数据
    frame = randi([0,1], frameLen, 1);
    syncPreSrc = frame(1 : syncLen);
    dataSrc = frame(syncLen+1 : frameLen-syncLen);
    syncPostSrc = frame(syncLen+dataLen+1 : frameLen);

%% 预编码
    frameCode = zeros(frameLen, 1);
    bitPre = 0;
    for n = 1 : frameLen/2
        if bitPre == 0 && frame(2*n-1) == 0
            bitPre = frame(2*n-1);
            frameCode(2*n-1) = 1;
        elseif bitPre == 1 && frame(2*n-1) == 0
            bitPre = frame(2*n-1);
            frameCode(2*n-1) = 0;
        elseif bitPre == 0 && frame(2*n-1) == 1
            bitPre = frame(2*n-1);
            frameCode(2*n-1) = 0;
        else
            bitPre = frame(2*n-1);
            frameCode(2*n-1) = 1;
        end

        if bitPre == 0 && frame(2*n) == 0
            bitPre = frame(2*n);
            frameCode(2*n) = 0;
        elseif bitPre == 1 && frame(2*n) == 0
            bitPre = frame(2*n);
            frameCode(2*n) = 1;
        elseif bitPre == 0 && frame(2*n) == 1
            bitPre = frame(2*n);
            frameCode(2*n) = 1;
        else
            bitPre = frame(2*n);
            frameCode(2*n) = 0;
        end
    end

%% 编码后数据
    syncPreCode = frameCode(1 : syncLen);
    dataCode = frameCode(syncLen+1 : frameLen-syncLen);
    syncPostCode = frameCode(syncLen+dataLen+1 : frameLen);
end
