function [voutestim,vref,g]=AGC(signal,N,time)
vref=0;
g=0;
Pow_sig_recvf=(1/N)*sum(abs(signal.^2));
voutestim=20*log10(Pow_sig_recvf);
disp(voutestim)
end

