function H = beamspace_channel(n,K,L)
% n: number of transmit beams (transmit antennas)
% K: number of users
% L: number of paths, L=1 for LoS, L>1 for multipath
lamada = 1; % wavelength
d = lamada/2; % distance between two antennas
H = zeros(n,K);
theta = zeros(L,K);
beta = zeros(L,K);
for k=1:K
    %%% complex gain
    beta(1:L,k) = (randn(1,L)+1i*randn(1,L))/sqrt(2);
    %%% DoA
    theta(1:L,k) = pi*rand(1,L) - pi/2;
    %%% channel for the k th user
    for j = 1:L
        H(:,k) = H(:,k) + beta(j,k)*array_respones(theta(j,k),n,d,lamada);
    end
end
% 0.5*sin(theta);
H = sqrt(n/L)*H;