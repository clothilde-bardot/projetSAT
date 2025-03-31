function [cn0, good1, good2, good_val1, good_val2, good_val11, good_val21, good_val12, good_val22, val] = beaulieu_v2(N,Tint,I)
    Pn=(abs(I(2:N))-abs(I(1:N-1))).^2;
    %Pn=(I(2:N)-I(1:N-1)).^2;
    Pd=1/2*(I(2:N).^2+I(1:N-1).^2);
    cn0=10*log10(1/(Tint)*(1/N*sum(Pn./Pd))^(-1));

    index = find(I>0);
    good1 = length(index); 
    % var_tt = Pn./Pd;
    % figure;
    % histfit(var_tt);
    % 
    % figure;
    % histfit(1./var_tt);

    Pn1 = ((I(2:N))-(I(1:N-1))).^2;
    good_val = (1/N*sum(Pn1./Pd));

    index1 = find(Pn ~= Pn1);
    index2 = find(Pn == Pn1);

    good2 = length(index1);
    good_val1 = mean(Pn(index1)./Pd(index1));
    good_val2 = mean(Pn(index2)./Pd(index2));

    good_val11 = mean(Pn(index1));
    good_val21 = mean(Pn(index2));

    good_val12 = mean(Pd(index1));
    good_val22 = mean(Pd(index2));

    val = 1/N*sum(Pn./Pd);


end