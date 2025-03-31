function [f_vec, H_filter, loss_pow_filter, Rh_real] = setting_RF_filter_from_BW_nom_passband_v2(BW, Fs, fIF, stepf_in, plot_active)

%génére une réponse fréquentielle pour un filtre RF et donne la réponse
%impulsionelle dans le domaine temporelle

limit = 100; %max frequency range in MHz
trans = 12.7; %stepf_in/1e6;% %transition width

%Generating vector freq
stepf = stepf_in/1e6; %step en MHz
f_vec = -Fs/2:stepf:Fs/2; %vecteur fréquentiel
len_vec = length(f_vec);
H_filter = -70*ones(1,len_vec); %initialisation réponse filtre à -70dB sur toutes les fréquences

if (trans == 0)
    slope = 0;
else
    slope = 70/trans;
end

%Definition in MHz
t1 = -limit:stepf:(-trans-BW/2-fIF); %avant la transition
t2 = (-trans-BW/2-fIF+stepf):stepf:(-BW/2-fIF-stepf); %début de la transition
t3 = -BW/2-fIF:stepf:BW/2-fIF; %passe bande
t4 = (BW/2-fIF+stepf):stepf:(BW/2-fIF+trans-stepf);%fin transition

t5 = (BW/2-fIF+trans):stepf:(-trans-BW/2+fIF);

t6 = (-trans-BW/2+fIF+stepf):stepf:(-BW/2+fIF-stepf);
t7 = -BW/2+fIF:stepf:BW/2+fIF;
t8 = (BW/2+fIF+stepf):stepf:(BW/2+fIF+trans-stepf);
t9 = (BW/2+fIF+trans):stepf:limit;

%réponse fréquentielle pour chaque plage fréquentielle
H_filter_temp = [-70+0*t1 -70+slope*(t2-(-BW/2-trans-fIF)) 0*t3 -slope*(t4-BW/2+fIF) -70+0*t5 -70+slope*(t6-(-BW/2-trans+fIF)) 0*t7 -slope*(t8-BW/2-fIF) -70+0*t9];
f_temp = [t1 t2 t3 t4 t5 t6 t7 t8 t9];
len_Htemp = length(f_temp);



f_vec2 = -(Fs+limit):stepf:(Fs+limit);
%ajuste H_filter_ramp pour etre sur toutes les fréquences
H_filter_temp2 = zeros(1, length(f_vec2));
%H_filter_temp2(1:len_Htemp) = H_filter_temp2(1:len_Htemp) + 10.^(H_filter_temp/20);
ind0 = find(abs(f_vec2)<1e-6);
ind01 = ind0 - floor(len_Htemp/2);
ind02 = ind0 + floor(len_Htemp/2);
H_filter_temp2(ind01:ind02) = H_filter_temp2(ind01:ind02) + 10.^(H_filter_temp/20);
%H_filter_temp2(end+1-len_Htemp:end) = H_filter_temp2(end+1-len_Htemp:end) + 10.^(H_filter_temp/20);
H_filter_temp2 = 20*log10(H_filter_temp2);
%H_filter_temp2(ind0) = (H_filter_temp2(ind0+1)+H_filter_temp2(ind0-1))/2;


indleft = find(abs(f_vec2 + Fs/2) < 1e-6);
indrigth = find(abs(f_vec2 - Fs/2) < 1e-6);

%troncature
H_filter = H_filter_temp2(indleft:indrigth);
f_vec = f_vec*1e6;


H_filter = H_filter(1:end-1);
f_vec = f_vec(1:end-1);

if (plot_active)
    figure;
    plot(f_vec/1e6, H_filter);
    title('fonction de transfert du filtre');
end
%calcule la puissance perdue (power loss) et la réponse impulsionelle en
%domaine temporel
H_lin = 10.^(H_filter/10);
H_lin_fft = fftshift(H_lin);
%Rh_time = ifft(H_lin_fft); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST
Rh_time = fftshift(ifft(H_lin_fft));
Rh_real = real(Rh_time);
Rh0 = max(Rh_real); % <-- I would say normalized power
Rh0_from_freq = sum(H_lin)/length(f_vec); % <-- I would say normalized power

loss_pow_filter = Rh0_from_freq;

if plot_active
    figure;
    plot(real(Rh_time));
    hold on;
    plot(imag(Rh_time));
    title('Réponse impulsionelle')
end

