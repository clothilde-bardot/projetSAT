function cn0 = snv(N,Tint,I,Q)

    Pd = (sum(I)/N).^2;
    Pdem = (sum(I.^2) + sum(Q.^2))/N;
    cn0=10*log10(1/(Tint)*(Pd/(Pdem-Pd)));

end
