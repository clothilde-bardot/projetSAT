function cn0 = NWPR(N,M, Tint,I,Q)

    nloop = N/M;

    NBP = zeros(1, nloop);
    WBP = zeros(1, nloop);
    NP = zeros(1, nloop);

    for i = 1:nloop

        NBP(i) = sum(I(1+(i-1)*M:i*M)).^2 + sum(Q(1+(i-1)*M:i*M)).^2;
        WBP(i) = sum(I(1+(i-1)*M:i*M).^2) + sum(Q(1+(i-1)*M:i*M).^2);
        NP(i) = NBP(i)/WBP(i);
    end

    NP_mean = mean(NP);

    cn0 = 10*log10(1/(Tint)*((NP_mean-1)/(M-NP_mean)));

end