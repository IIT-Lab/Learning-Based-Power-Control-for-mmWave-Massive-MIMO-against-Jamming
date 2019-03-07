function C3 = main_MU(P_B,P_J,N_T)
SNR_dB=5;  SNR_linear=10.^(SNR_dB/10.);
N_iter=1; 
%N_T = 64;
N_J = 16;

K = 16;
Nt_RF= K;
L = 3; % number of rays
lamada = 1; % wavelegenth;


for i_snr=1:length(SNR_linear)
    
    SNR=SNR_linear(i_snr);
    temp=0;temp1=0;temp2=0;temp3=0;temp4=0;temp5=0;
    for i=1:N_iter
        H_T = beamspace_channel(N_T,K,L).';
        H_J = beamspace_channel(N_J,K,L).';

        %%% ACE-based hybrid precoding
        F3 = ACE_sub(H_T,SNR,Nt_RF,1,K);
        F2 = ACE_sub(H_J,SNR,Nt_RF,1,K);
        beta3 = sqrt(K/trace(F3*F3'));
        beta2 = sqrt(K/trace(F2*F2'));
        H_eq3=H_T*F3;
        G_eq3=H_J*F2;
        for k=1:K
            sum_inf = sum(abs(H_eq3(:,k)).^2)-abs(H_eq3(k,k))^2;
            sumG = sum(abs(G_eq3(:,k)).^2);
            temp3 = temp3+log2(1+P_B*abs(H_eq3(k,k))^2/(P_B*sum_inf+P_J*sumG+K/(SNR*beta3^2)));
        end      
        
        %%% AS-based hybrid precoding
%         F4 = AS(H,Nt_RF);
%         beta4 = sqrt(K/trace(F4*F4'));
%         H_eq4=H*F4;
%         for k=1:K
%             sum_inf = sum(abs(H_eq4(:,k)).^2)-abs(H_eq4(k,k))^2;
%             temp4 = temp4 + log2(1+abs(H_eq4(k,k))^2/(sum_inf+K/(SNR*beta4^2)));
%         end    
        
        
    end
%     C(i_snr)= real(temp/N_iter);
%     C1(i_snr)= real(temp1/N_iter);
%     C2(i_snr)= real(temp2/N_iter);     
    C3(i_snr)= real(temp3/N_iter);
%     C4(i_snr)= real(temp4/N_iter);
end
end
