function F = two_stage(H,Nt_RF,bit)
K = size(H,1);
N = size(H,2);
A = zeros(N,Nt_RF);
N_phase = 2^bit;
phase_set = (2*pi/N_phase)*[0:N_phase-1] - pi;
for k = 1:K
    phase = angle(H(k,:)');
    for n = 1:N
        [~,order] = min(abs(phase(n) - phase_set));
        Q_phase(n,1) = phase_set(order);
    end
    A(:,k) = exp(1i*Q_phase);
end
H_eq = H*A;
D = H_eq'*inv(H_eq*H_eq');
F = A*D;