function A = sub_matrix(B,M,Nt_RF)
N = M*Nt_RF;
C = reshape(B,M,Nt_RF);
A = zeros(N,Nt_RF);
for i = 1 : Nt_RF
    A(M*(i-1)+[1:M],i) = C(:,i);
end