function [cn0, M2, M4, M4v3, M41, M42, M43] = MM_v2(N,Tint,I,Q)
    M2=1/N*sum(I.^2 + Q.^2);
    M4=1/N*sum((I.^2 + Q.^2).^2);
    M41 = 1/N*sum(I.^4);
    M42 = 1/N*sum(Q.^4);
    M43 = 2/N*sum((I.^2).*(Q.^2));
    M4v3 = M41 + M42 + M43;
    Pd=sqrt(2*M2^2-M4);
    cn0=10*log10(1/(2*Tint)*2*Pd/(M2-Pd));
end