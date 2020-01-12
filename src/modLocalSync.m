function localSyncModulated = modLocalSync(localSync)
% localSync : ����ͬ������
% localSyncModulated : ���ƺ�ı���ͬ������
    len = length(localSync);
    localSyncModulated = zeros(len-1, 1); 
    for index = 1 : 2 : len-1
        localSyncModulated(index) = localSync(index+1) + 1j * localSync(index);
    end
    for index = 2 : 2 : len-1
        localSyncModulated(index) = localSync(index) + 1j * localSync(index+1);
    end  
end
