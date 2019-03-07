function F = AS(H,Nt_RF)
K = size(H,1);
N = size(H,2);
deta = 10^-1;
% H£ºbeamspace channel matrix (sparse)
% N; number of retained beams 
% n: total number of beams (transmit antennas)
H_abs = abs(H);
s_set = 1:N;
T_set = [];
%%% proposed search 
for i = 1:Nt_RF
    A = H(:,T_set);
    B = inv(A*A' + deta*eye(K,K));
    target = [];
    for j = 1:length(s_set)
        temp = H(:,s_set(j));
        if ~isempty(A)
           target(j) = (temp'*B^2*temp)/(1+temp'*B*temp);
        else
           target(j) = (temp'*temp)/(1+temp'*temp);
        end
    end
    [~,jj] = max(target);
    T_set = unique([s_set(jj) T_set]);
    s_set(jj) = [];
end
S = zeros(N,Nt_RF);
for i = 1:Nt_RF
    S(T_set(i),i) = 1;
end
% H_eff1 = H*S;
H_eff = H(:,T_set);
D = H_eff'*inv(H_eff*H_eff');
F = S*D;
