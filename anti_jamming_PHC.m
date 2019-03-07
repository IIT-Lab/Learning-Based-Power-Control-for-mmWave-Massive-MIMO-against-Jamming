%%%%%%%%%%%%%%%%%%anti_jamming
function anti_jamming_PHC()
clc;
clear all;



%  load Q_tr.mat
Ct = 1.5;    %PHC��q�õĲ���Ct = 0.1;Cj = 0.8;
Cj = 5;
num = 0;
N_antenna = [48,64,96,128,160,256];
for N_liner = 1:length(N_antenna)
    datestr(now,0)
    N = N_antenna(N_liner)
    num = num + 1;
% H_eq3 = H_BSinitial(SNR_dB);
% G_eq3 = H_BSinitial(SNR_dB) +0.1;% K = [0.639785973030503,0.302280985837198,0.498598218246061,0.325251136513776;0.812070237688966,0.130808261405803,0.804978312373026,0.892138207548114;0.156571481194958,0.349831708722582,0.816503254382907,0.625626614086431;0.825102566581947,0.836965539227896,0.0707019263805208,0.372754374902034];


transmitter_R_space = 0:40;
transmitter_action_space = [1,2,3,4,5,6,7,8,9,10];
transmitter_power_space = [1,2,3,4,5,6,7,8,9,10];
jammer_action_space = [1,2,3,4,5,6,7,8,9,10];
jammer_state_space = [1,2,3,4,5,6,7,8,9,10];

alpha1 = 0.5;
delta1 = 0.85;
alpha2 = 0.5;
delta2 = 0.5;

num_episode = 30;
time_step = 5000;
epsilon1 = 1;   %PHC
epsilon2 = 0.1;

% sigma = 0;
% sigma_w = 0.03;   %PHCƫ��
% sigma_s = 0.05;
w = 0.09;
transmitter_R_num = length(transmitter_R_space);
transmitter_action_num = length(transmitter_action_space);
transmitter_power_num = length(transmitter_power_space);
jammer_action_num = length(jammer_action_space);
jammer_state_num = length(jammer_state_space);



for episode = 1:num_episode

    episode
    transmitter_power_selection = transmitter_power_space(randi(transmitter_power_num));
     transmitter_R_selection = randi(transmitter_R_num);
    jammer_state_selection = jammer_state_space(randi(jammer_state_num));
    
%     Q_transmitter = Q_; 
    Q_transmitter = zeros(transmitter_R_num,transmitter_power_num,transmitter_action_num);
    Q_jammer = zeros(jammer_state_num,jammer_action_num);
    
    posibility_matrix = ones(transmitter_R_num,transmitter_power_num,transmitter_action_num)/10;
%     ones(transmitter_R_num,transmitter_power_num,transmitter_action_num)/10;   %PHC strategy  ���ʶ���ѡ�����?
%     average_posibility_matrix = ones(transmitter_state_num,transmitter_action_num)/10;
%     Cs = zeros(1,length(transmitter_state_space));
    
    for step = 1:time_step
   
%          transmitter_state_record(episode,step) = transmitter_R_selection;
        %%transmitter
        matrix = reshape(posibility_matrix(transmitter_R_selection,transmitter_power_selection,:),1,10);
        transmitter_action_selection = ...
            randsrc(1,1,[transmitter_action_space;matrix]);
        
%         %greedy
%         possibility1 = rand();
%         if possibility1 < (1 - epsilon1)
%             Q_transmitter_max = max(Q_transmitter(transmitter_state_selection,:));
%             index1 = find(Q_transmitter(transmitter_state_selection,:) == Q_transmitter_max);
%             transmitter_next_action_selection = transmitter_action_space(index1(randi(length(index1)))); 
%         else
%             transmitter_next_action_selection = transmitter_action_space(randi(transmitter_action_num));  
%         end
        
        %jammer
        possibility2 = rand();
        if possibility2 < (1 - epsilon2)
            Q_jammer_max = max(Q_jammer(jammer_state_selection,:));
            index2 = find(Q_jammer(jammer_state_selection,:) == Q_jammer_max);
            jammer_next_action_selection = jammer_action_space(index2(randi(length(index2))));  
        else
            jammer_next_action_selection = jammer_action_space(randi(jammer_action_num));   
        end
        

        
        
        %transmitter_utility
        transmitter_rate = main_MU(transmitter_action_selection,jammer_next_action_selection,N);%rateFunc(G_eq3,H_eq3,transmitter_action_selection,jammer_next_action_selection,SNR_dB); 

        transmitter_utility = transmitter_rate - Ct*transmitter_action_selection;

        %jammer_utility  
%         jammer_utility = - transmitter_utility;
        jammer_rate = -main_MU(transmitter_action_selection,jammer_next_action_selection,N);

        jammer_utility = jammer_rate - Cj*jammer_next_action_selection;%-transmitter_utility;
            
        
        
        %update Q_table        
        Q_transmitter(transmitter_R_selection,transmitter_power_selection,transmitter_action_selection)=...
            (1-alpha1)*Q_transmitter(transmitter_R_selection,transmitter_power_selection,transmitter_action_selection)+...
            alpha1*(transmitter_utility+delta1*max(Q_transmitter(transmitter_R_selection,transmitter_power_selection,:)));
        
        Q_jammer(jammer_state_selection,jammer_next_action_selection)=...
            (1-alpha2)*Q_jammer(jammer_state_selection,jammer_next_action_selection)+...
            alpha2*(jammer_utility+delta2*max(Q_jammer(jammer_state_selection,:)));   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %����ƽ����ԵĹ���?%         Cs(transmitter_state_selection) =  Cs(transmitter_state_selection) + 1;
%         average_posibility_matrix(transmitter_state_selection,:) = average_posibility_matrix(transmitter_state_selection,:) + (posibility_matrix(transmitter_state_selection,:) - average_posibility_matrix(transmitter_state_selection,:))/Cs(transmitter_state_selection);
%         if Q_transmitter(transmitter_state_selection,:) * (posibility_matrix(transmitter_state_selection,:)') > Q_transmitter(transmitter_state_selection,:) * (average_posibility_matrix(transmitter_state_selection,:)')
%             w = sigma_w;
%         else
%             w = sigma_s;
%         end
        %����probablity_matrix ��
        [~,transmitter_next_action_selection_temp] = max(Q_transmitter(transmitter_R_selection,transmitter_power_selection,:));
        posibility_matrix(transmitter_R_selection,transmitter_power_selection,transmitter_next_action_selection_temp) = ...
            posibility_matrix(transmitter_R_selection,transmitter_power_selection,transmitter_next_action_selection_temp) + w;
        if posibility_matrix(transmitter_R_selection,transmitter_power_selection,transmitter_next_action_selection_temp) >= 1
            posibility_matrix(transmitter_R_selection,transmitter_power_selection,transmitter_next_action_selection_temp) = epsilon1;
            for action_temp = 1:length(transmitter_action_space)
                if action_temp ~= transmitter_next_action_selection_temp
                    posibility_matrix(transmitter_R_selection,transmitter_power_selection,action_temp) = (1 - epsilon1)/(length(transmitter_action_space) - 1);
                end
            end
        else
            sum_temp = 0;
            for action_temp = 1:length(transmitter_action_space)
                if action_temp ~= transmitter_next_action_selection_temp
                    posibility_matrix(transmitter_R_selection,transmitter_power_selection,action_temp) = ...
                        posibility_matrix(transmitter_R_selection,transmitter_power_selection,action_temp) - w/(length(transmitter_action_space) - 1);
                    if posibility_matrix(transmitter_R_selection,transmitter_power_selection,action_temp) < 0
                        sum_temp = sum_temp + posibility_matrix(transmitter_R_selection,transmitter_power_selection,action_temp);
                        posibility_matrix(transmitter_R_selection,transmitter_power_selection,action_temp) = 0;
                    end
                end
            end
            posibility_matrix(transmitter_R_selection,transmitter_power_selection,transmitter_next_action_selection_temp) = posibility_matrix(transmitter_R_selection,transmitter_power_selection,transmitter_next_action_selection_temp) + sum_temp;
        end
                        
        %update
        transmitter_power_selection = jammer_next_action_selection;
        [~,transmitter_R_selection] = find(transmitter_R_space == floor(transmitter_rate));
        jammer_state_selection = transmitter_action_selection;
        
        
        %record
        transmitter_utility_record(num,episode,step) = transmitter_utility;
        jammer_utility_record(num,episode,step) = jammer_utility;
         transmitter_rate_record(num,episode,step) = transmitter_rate;
%          SINR_record(num,episode,step) = 2^transmitter_rate - 1;
%         transmitter_state_record(episode,step) = transmitter_state_selection;
    
    end
    
end
end
save('data_PHC.mat');

%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%
% transmitter_state_average = mean(transmitter_state_record);
SINR_record = (2.^(transmitter_rate_record/16)-1)*16;
for index = 1:6
    sum_rate_PHC(index) = sum(mean(reshape(transmitter_rate_record(index,:,:),num_episode,time_step)))/time_step;
    sum_sinr_PHC(index) = sum(mean(reshape(SINR_record(index,:,:),num_episode,time_step)))/time_step;
    sum_utility_PHC(index) = sum(mean(reshape(transmitter_utility_record(index,:,:),num_episode,time_step)))/time_step;
end

transmitter_utility_recordd = reshape(transmitter_utility_record(3,:,:),num_episode,time_step);
% jammer_utility_record = reshape(jammer_utility_record(7,:,:),200,10000);
transmitter_rate_recordd = reshape(transmitter_rate_record(3,:,:),num_episode,time_step);
 SINR_recordd = reshape(SINR_record(3,:,:),num_episode,time_step);

 transmitter_utility_average = mean(transmitter_utility_recordd);
% jammer_utility_average = mean(jammer_utility_record);
 transmitter_rate_average = mean(transmitter_rate_recordd);
 SINR_average = mean(SINR_recordd);
% 
 transmitter_utility_averag = smooth(transmitter_utility_average,30);
% jammer_utility_averag = smooth(jammer_utility_average,30,'loess');
 transmitter_rate_averag = smooth(transmitter_rate_average,30);
 SINR_averag = smooth(SINR_average,30);
 
figure(1);
plot(transmitter_utility_averag,'k','linewidth',1);
xlabel('Time slot');
ylabel('Utility');
legend('PHC','Q-learning','RCRA');
grid on;
hold on;

% figure(2);
% plot(jammer_utility_averag);
% xlabel('time slot');
% ylabel('jammer utility');
% hold on;

figure(2);
plot(transmitter_rate_averag,'k','linewidth',1);
xlabel('Time slot');
ylabel('Sum data rate');
legend('PHC','Q-learning','RCRA');
grid on;
hold on;

figure(3);
plot(SINR_averag,'k','linewidth',1);
xlabel('Time slot');
ylabel('Sum SINR');
legend('PHC','Q-learning','RCRA');
grid on;
hold on;

% w = 16:16:96;
figure(4);
plot(N_antenna,sum_rate_PHC,'-ks','linewidth',1);
xlabel('The number of the antennas');
ylabel('Sum date rate');
legend('PHC','Q-learning','RCRA');
grid on;
hold on;

figure(5);
plot(N_antenna,sum_utility_PHC,'-ks','linewidth',1);
xlabel('The number of the antennas');
ylabel('Average utility');
legend('PHC','Q-learning','RCRA');
grid on;
hold on;

figure(6);
plot(N_antenna,sum_sinr_PHC,'-ks','linewidth',1);
xlabel('The number of the antennas');
ylabel('Average sum SINR');
legend('PHC','Q-learning','RCRA');
grid on;
hold on;


end



























