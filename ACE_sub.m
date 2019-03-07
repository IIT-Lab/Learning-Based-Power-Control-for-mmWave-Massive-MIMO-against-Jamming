function [F,A_opt,D_opt]=ACE_sub(H,SNR,Nt_RF,N_iter,K)
U = size(H,1);
N = size(H,2);
M = N/Nt_RF;
% K = 200;
% T = N * Nt_RF;
rho = 0.2;
% p(1,:) = (1/2)*ones(1,T);
% rho = 0.2;
p(1,:) = (1/2)*ones(1,N);
R = ceil((1-rho)*K);
S = zeros(1,K);
for iter = 1:N_iter
    for j = 1:N
        X(:,j) = randsrc(K,1,[1,-1; p(iter,j), 1 -  p(iter,j)]);
    end
    for k = 1:K
        A = sub_matrix(X(k,:),M,Nt_RF);
        
%         A = reshape(X(k,:),N,Nt_RF);
%         if rank(A) < Nt_RF
%           continue
%         end
        H_eff = H*A;
        D = H_eff'*inv(H_eff*H_eff');
        F = A*D;
        H_eq = H*F;
        beta = sqrt(K/trace(F*F'));
        temp = 0;
        for u=1:U
            sum_inf = sum(abs(H_eq(:,u)).^2)-abs(H_eq(u,u))^2;
            temp = temp + log2(1+abs(H_eq(u,u))^2/(sum_inf+K/(SNR*beta^2)));
        end
        S(k) = temp;
    end
    [S1,order] = sort(S,'ascend');
    r(iter) = S1(R);
    candidate = find(S>=r(iter));
    T = sum(S(candidate))/length(candidate);
    den = length(candidate);
    for j = 1:N
        num = 0;
        for m = 1:den
            if X(candidate(m),j) == 1
                num = num + S(candidate(m))/T;
            end
        end
        p(iter+1,j) = num/den;
        if p(iter+1,j) > 1
            p(iter+1,j) = 1;
        end
    end
end
A_opt = sub_matrix(X(order(end),:),M,Nt_RF);
% A_opt = reshape(X(order(end),:),N,Nt_RF);
H_eff = H*A_opt;
D_opt = H_eff'*inv(H_eff*H_eff');
F = A_opt*D_opt;