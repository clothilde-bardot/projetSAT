
%INITIALISATION
clear all;
close all;
clc;

%PARAMÈTRES GÉNÉRAUX
%rng(1001);
nb_loops = 10;
val_A = 0.5; %à voir avec var_sigma %%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%%%%
%val_A = 0.58;
nb_bits = 3;


CN0 = 40;
var_sigma = 1;
Ti = 1e-3;
Fs = 100e6;   %fréquence d'échantillonage
PRN = 10;   %code pseudo-aléatoire GPS, indice ds la liste
fIF = 25e6; %fréquence intermédiaire (cf filtre bandwith)

BW = 20; %in MHz %bande passante
stepf = 1e3; %in Hz



%PARAMETRES AJOUTES (AGC)
L = 2^(nb_bits-1) - 1;    % Offset de quantification (3 pour 3 bits)
B = 2^nb_bits - 1;        % Plage max (7 pour 3 bits)
target_level = 0.6 * L;       % Cible: 80% de la plage de quantification
alpha = 0.05;                % Coefficient de lissage (gain de la boucle)
window_size = 100;           % Taille de la fenêtre pour le calcul RMS
val_Aq = val_A;              % Initialisation du gain AGC (valeur de départ = val_A)




%PARAMÈTRES DU SIGNAL GNSS
fc = 10.23e6;
Tc = 1/fc;
Ts = 1/Fs;
N = Ti/Ts;
N1s = 1/Ts; % échantillons sur 1 sec pour estimation C/No
M = 1/Ti; % nmb de sorties du corrélateurs sur 1s

%GÉNÉRATION DU CODE PRN
[code_L5I, code_L5Q] = GPS_L5_code_prim(PRN); %cf code : donne voies I et Q modulation BPSK

%FILTRAGE : => donne la réponse fréquentielle du filtre, H_filter = gain du
%filtre en dB
% filter generation
[f_vec, H_filter, loss_pow_filter, Rh_real] = setting_RF_filter_from_BW_nom_passband_v2(BW, Fs/1e6, fIF/1e6, stepf, 1);
H_lin_sqrt = 10.^(H_filter/20);

%paramétres de quantification
L = 2^(nb_bits-1) - 1; %nb de bits de quantif
B = 2^nb_bits - 1; %niveaux de quantification

%GÉNÉRATION DU SIGNAL GNSS SIMULÉ
%génére un signal modulé
time = (0:N1s-1)*Ts;
time_tc = time/Tc; %normalisation du temps par Tc (Temps d'un chip)
len_time = length(time);
code_I = create_code_samples(code_L5I, time_tc); %cf code (échantillone le code PRN sur le vecteur de temps voulu)
carrier_I = cos(2*pi*fIF*time); %signal I en phase
carrier_Q = sin(2*pi*fIF*time); %signal Q en quadrature
Amp = sqrt(4*10^(CN0/10)*Ts); %amplitude du signal en fct de CNo
sig_useful = Amp*code_I.*carrier_I; %signal modulé (le code i module l'amplitude de la porteuse en phase I) (cf état de l'art équation(3))

%estimation des CNo

CNO_est = zeros(1, nb_loops);
CNO_estf = zeros(1, nb_loops);
CNO_estfq = zeros(1, nb_loops);
CNO_estq = zeros(1, nb_loops);

CNO_NWPR = zeros(1, nb_loops);
CNO_NWPRf = zeros(1, nb_loops);
CNO_NWPRq = zeros(1, nb_loops);
CNO_NWPRfq = zeros(1, nb_loops);
CNO_NWPR2 = zeros(1, nb_loops);
CNO_NWPRf2 = zeros(1, nb_loops);

CNO_SNV = zeros(1, nb_loops);
CNO_SNVf = zeros(1, nb_loops);
CNO_SNVq = zeros(1, nb_loops);
CNO_SNVfq = zeros(1, nb_loops);

CNO_MM = zeros(1, nb_loops);
CNO_MMf = zeros(1, nb_loops);
CNO_MMq = zeros(1, nb_loops);
CNO_MMfq = zeros(1, nb_loops);

CNO_B = zeros(1, nb_loops);
CNO_Bf = zeros(1, nb_loops);
CNO_Bq = zeros(1, nb_loops);
CNO_Bfq = zeros(1, nb_loops);

psd_tot = zeros(1, N);

%simulation de la réception, filtrage, quantification et corrélation entre
%le signal recu et le code PNR

%BOUCLE DE SIMULATION

for i = 1:nb_loops %nb d'itération
    disp(['Num Loops = ' num2str(i) '/' num2str(nb_loops)])
    
    %INITIALISATION DES VECTEURS : matrices pour stocker les corrélation du
    %signal etdu bruit

    %vecteur de corrélation du signal recu
    corr_I = zeros(1, M); 
    corr_Q = zeros(1, M);
    %vecteur de corrélation du signal utile
    corr_use_I = zeros(1, M);
    corr_use_Q = zeros(1, M);
    %bruit
    corr_noise_I = zeros(1, M);
    corr_noise_Q = zeros(1, M);
    %signal filtré
    corr_If = zeros(1, M);
    corr_Qf = zeros(1, M);
    %signal utile filtré
    corr_use_If = zeros(1, M);
    corr_use_Qf = zeros(1, M);
    %bruit filtré
    corr_noise_If = zeros(1, M);
    corr_noise_Qf = zeros(1, M);
    %signal quantifié
    corr_Iq = zeros(1, M);
    corr_Qq = zeros(1, M);
    %signal quantifié et filtré
    corr_Ifq = zeros(1, M);
    corr_Qfq = zeros(1, M);

    rms_output_history = zeros(1, M); % Stockage des valeurs RMS pour chaque itération
    gain_history = zeros(1, M);       % Historique du gain val_Aq
    error_history = zeros(1, M);      % Historique de l'erreur
    output_voltage_history = zeros(1, M); % Historique de la tension de sortie

    
    for n = 1:M
        
        %SIMULATION DU BRUIT (WGN)
        noise_I   = randn(1, len_time/M); %génération du bruit (bruit gaussien blanc)
        %noise_Q   = randn(1, len_time);
     
        noise_If   = real(ifft(fftshift(fftshift(fft(noise_I)).*H_lin_sqrt))); %filtrage du bruit par un passe bande
        pow_noisef = var(noise_If);

        indini = 1 + (n-1)*N;
        indend = N + (n-1)*N;

        time_local = time(indini:indend);
        
        %SIMULATION DU SIGNAL RECU
        sig_use_recv = sig_useful(indini:indend); %signal utile sur une fenêtre temporelle considérée avec sig_useful = Amp*code_I.*carrier_I;
        sig_use_recvf = real(ifft(fftshift(fftshift(fft(sig_use_recv)).*H_lin_sqrt)));%filtrage de ce signal

        sig_recv = sig_use_recv + noise_I; %signal recu plus le bruit
        sig_recvf = sig_use_recvf + noise_If; % signal recu plus le bruit mais filtré
        
        %QUANTIFICATION DU SIGNAL
        %Quantizer - No filter
        %val_Aq = val_A;
        %temp_vec_q = max(0, ceil(val_Aq*sig_recv+L)); %amplitude =>niveau de quantification
        %sig_recvq = min(B, -B + 2*temp_vec_q); %quantification du signal recu


        % Quantizer - No filter AJOUT AGC
        % --- AGC DYNAMIQUE ---
        % 1. Application du gain
        sig_recv_scaled = val_Aq * sig_recv;

        % 2. Quantification
        sig_recvq = min(B, -B + 2 * max(0, ceil(sig_recv_scaled + L)));

        % 3. Mise à jour du gain (si on a assez d'échantillons)
        if n > window_size
            % Mesure RMS du signal quantifié
            rms_output = sqrt(mean(sig_recvq(n-window_size:n).^2));
            rms_output_history(n) = rms_output; % Sauvegarde de la valeur RMS
    
            % Correction du gain (Loop Filter)
            error = target_level - rms_output;
            val_Aq = val_Aq + alpha * error;
    
            % Limites de sécurité
            % val_Aq = max(0.1, min(2.0, val_Aq));

        % Enregistrement des valeurs à chaque itération
        gain_history(n) = val_Aq;
        error_history(n) = target_level - rms_output;
        output_voltage_history(n) = mean(sig_recvq(n-min(window_size,n-1):n)); % Moyenne glissante



        end

        %Quantizer - Filter
        val_Afq = val_A/sqrt(pow_noisef);%prend en compte la puissance du bruit filtré
        temp_vec_fq = max(0, ceil(val_Afq*sig_recvf+L));
        sig_recvfq = min(B, -B + 2*temp_vec_fq);

        %psd_tot = psd_tot + abs(fftshift(fft(sig_recvfq))).^2/N;
        %psd_tot = psd_tot + abs(fftshift(fft(sig_use_recv))).^2/N;

        carrier_local_I = carrier_I(indini:indend); %porteuse en phase sur l'intervalle
        carrier_local_Q = carrier_Q(indini:indend); %porteuse en quadrature sur l'intervalle
        code_local = code_I(indini:indend); %extrait le code PRN sur une fenêtre

        %CORRELATION ENTRE LE SIGNAL ET LE CODE PRN
        corr_I(n) = real(sum(sig_recv.*carrier_local_I.*code_local))/N; %corr entre code PRN et I
        corr_Q(n) = real(sum(sig_recv.*carrier_local_Q.*code_local))/N; %corr entre code PRN et Q
        corr_If(n) = real(sum(sig_recvf.*carrier_local_I.*code_local))/N; %mm chose mais filtré
        corr_Qf(n) = real(sum(sig_recvf.*carrier_local_Q.*code_local))/N;
        corr_Iq(n) = real(sum(sig_recvq.*carrier_local_I.*code_local))/N; %mm chose mais quantifié
        corr_Qq(n) = real(sum(sig_recvq.*carrier_local_Q.*code_local))/N;
        corr_Ifq(n) = real(sum(sig_recvfq.*carrier_local_I.*code_local))/N; %filtré+quantifié
        corr_Qfq(n) = real(sum(sig_recvfq.*carrier_local_Q.*code_local))/N;

        corr_use_I(n) = real(sum(sig_use_recv.*carrier_local_I.*code_local))/N; %pour le signal utile
        corr_use_Q(n) = real(sum(sig_use_recv.*carrier_local_Q.*code_local))/N;
        corr_use_If(n) = real(sum(sig_use_recvf.*carrier_local_I.*code_local))/N;
        corr_use_Qf(n) = real(sum(sig_use_recvf.*carrier_local_Q.*code_local))/N;
        
        corr_noise_I(n) = real(sum(noise_I.*carrier_local_I.*code_local))/N; % pour le bruit
        corr_noise_If(n) = real(sum(noise_If.*carrier_local_I.*code_local))/N;
        corr_noise_Q(n) = real(sum(noise_I.*carrier_local_Q.*code_local))/N;
        corr_noise_Qf(n) = real(sum(noise_If.*carrier_local_Q.*code_local))/N;


        % Afficher le gain AGC (optionnel)
        if mod(n, 100) == 0
            disp(['n = ', num2str(n), ' | val_Aq = ', num2str(val_Aq)]);
        end


    end

    %psd_tot = psd_tot/M;

    %ESTIMATION DU C/NO

    CN0_est(i) = mean(corr_use_I.^2)/mean(corr_noise_I.^2);
    CN0_est(i) = 10*log10(CN0_est(i)/2/Ti);

    CN0_estf(i) = mean(corr_use_If.^2)/mean(corr_noise_If.^2);
    CN0_estf(i) = 10*log10(CN0_estf(i)/2/Ti);

    CN0_estq(i) = mean(corr_Iq.^2)/var(corr_Iq);
    CN0_estq(i) = 10*log10(CN0_estq(i)/2/Ti);

    CN0_estfq(i) = mean(corr_Ifq.^2)/var(corr_Ifq);
    CN0_estfq(i) = 10*log10(CN0_estfq(i)/2/Ti);

    CN0_NWPR(i) = NWPR(M, 20, Ti, corr_I, corr_Q);
    CN0_NWPR2(i) = NWPR(1000, 1000, Ti, corr_I, corr_Q);
    CN0_NWPRf(i) = NWPR(M, 20, Ti, corr_If, corr_Qf);
    CN0_NWPRf2(i) = NWPR(1000, 1000, Ti, corr_If, corr_Qf);
    CN0_NWPRq(i) = NWPR(1000, 1000, Ti, corr_Iq, corr_Qq);
    CN0_NWPRfq(i) = NWPR(1000, 1000, Ti, corr_Ifq, corr_Qfq);

    CN0_SNV(i) = snv(M, Ti, corr_I, corr_Q);
    CN0_SNVf(i) = snv(M, Ti, corr_If, corr_Qf);
    CN0_SNVq(i) = snv(M, Ti, corr_Iq, corr_Qq);
    CN0_SNVfq(i) = snv(M, Ti, corr_Ifq, corr_Qfq);

    CN0_MM(i) = MM_v2(M, Ti, corr_I, corr_Q);
    CN0_MMf(i) = MM_v2(M, Ti, corr_If, corr_Qf);
    CN0_MMq(i) = MM_v2(M, Ti, corr_Iq, corr_Qq);
    CN0_MMfq(i) = MM_v2(M, Ti, corr_Ifq, corr_Qfq);

    CN0_B(i) = beaulieu_v2(M, Ti, corr_I);
    CN0_Bf(i) = beaulieu_v2(M, Ti, corr_If);
    CN0_Bq(i) = beaulieu_v2(M, Ti, corr_Iq);
    CN0_Bfq(i) = beaulieu_v2(M, Ti, corr_Ifq);


end

disp(['val_Aq = ', num2str(val_Aq)]);

figure; plot(rms_output_history); hold on; 
yline(target_level, 'r--', 'Target'); 
title('Amplitude Efficace du Signal Quantifié');

% Création d'une figure avec 3 sous-graphiques
figure;

% Courbe 1 : Évolution du gain
subplot(3,1,1);
plot(gain_history, 'LineWidth', 1.5);
title('Évolution du gain AGC (val\_Aq)');
xlabel('Itérations');
ylabel('Gain');
grid on;

% Courbe 2 : Évolution de l'erreur
subplot(3,1,2);
plot(error_history, 'LineWidth', 1.5);
hold on;
yline(0, 'r--', 'LineWidth', 1); % Ligne zéro pour référence
title('Erreur (Target - RMS Output)');
xlabel('Itérations');
ylabel('Erreur');
grid on;

% Courbe 3 : Tension de sortie après AGC
subplot(3,1,3);
plot(output_voltage_history, 'LineWidth', 1.5);
yline(target_level, 'r--', 'LineWidth', 1); % Ligne de référence
title('Tension de sortie après AGC');
xlabel('Itérations');
ylabel('Amplitude');
grid on;

% Ajustement de l'espacement
sgtitle('Analyse des Performances de l''AGC');


CN0_est_total = mean(CN0_est);
CN0_est_totalf = mean(CN0_estf);
CN0_est_totalq = mean(CN0_estq);
CN0_est_totalfq = mean(CN0_estfq);
disp(['CN0_est = ' num2str(CN0_est_total) ' -- CN0_estf = ' num2str(CN0_est_totalf) ' -- CN0_estq = ' num2str(CN0_est_totalq) ' -- CN0_estfq = ' num2str(CN0_est_totalfq)]);

CN0_NWPR_total = mean(CN0_NWPR);
CN0_NWPR_total2 = mean(CN0_NWPR);
CN0_NWPR_totalf = mean(CN0_NWPRf);
CN0_NWPR_totalf2 = mean(CN0_NWPRf);
CN0_NWPR_totalq = mean(CN0_NWPRq);
CN0_NWPR_totalfq = mean(CN0_NWPRfq);
disp(['CN0_NWPR = ' num2str(CN0_NWPR_total) ' -- CN0_NWPRf = ' num2str(CN0_NWPR_totalf) ' -- CN0_NWPRq = ' num2str(CN0_NWPR_totalq) ' -- CN0_NWPRfq = ' num2str(CN0_NWPR_totalfq)]);
disp(['diffCN0_NWPRf = ' num2str(CN0_NWPR_total-CN0_NWPR_totalf) ' -- diffCN0_NWPRq = ' num2str(CN0_NWPR_total-CN0_NWPR_totalq) ' -- diffCN0_NWPRfq = ' num2str(CN0_NWPR_total-CN0_NWPR_totalfq)]);

CN0_SNV_total = mean(CN0_SNV);
CN0_SNV_totalf = mean(CN0_SNVf);
CN0_SNV_totalq = mean(CN0_SNVq);
CN0_SNV_totalfq = mean(CN0_SNVfq);
%disp(['CN0_SNV = ' num2str(CN0_SNV_total) ' -- CN0_SNVf = ' num2str(CN0_SNV_totalf) ' -- CN0_SNVq = ' num2str(CN0_SNV_totalq) ' -- CN0_SNVfq = ' num2str(CN0_SNV_totalfq)]);
disp(['diffCN0_SNVf = ' num2str(CN0_SNV_total-CN0_SNV_totalf) ' -- diffCN0_SNVq = ' num2str(CN0_SNV_total-CN0_SNV_totalq) ' -- diffCN0_SNVfq = ' num2str(CN0_SNV_total-CN0_SNV_totalfq)]);

CN0_MM_total = mean(CN0_MM);
CN0_MM_totalf = mean(CN0_MMf);
CN0_MM_totalq = mean(CN0_MMq);
CN0_MM_totalfq = mean(CN0_MMfq);
%disp(['CN0_MM = ' num2str(CN0_MM_total) ' -- CN0_MMf = ' num2str(CN0_MM_totalf) ' -- CN0_MMq = ' num2str(CN0_MM_totalq) ' -- CN0_MMfq = ' num2str(CN0_MM_totalfq)]);
disp(['diffCN0_MMf = ' num2str(CN0_MM_total-CN0_MM_totalf) ' -- diffCN0_MMq = ' num2str(CN0_MM_total-CN0_MM_totalq) ' -- diffCN0_MMfq = ' num2str(CN0_MM_total-CN0_MM_totalfq)]);

CN0_B_total = mean(CN0_B);
CN0_B_totalf = mean(CN0_Bf);
CN0_B_totalq = mean(CN0_Bq);
CN0_B_totalfq = mean(CN0_Bfq);
%disp(['CN0_B = ' num2str(CN0_B_total) ' -- CN0_Bf = ' num2str(CN0_B_totalf) ' -- CN0_Bq = ' num2str(CN0_B_totalq) ' -- CN0_Bfq = ' num2str(CN0_B_totalfq)]);
disp(['diffCN0_Bf = ' num2str(CN0_B_total-CN0_B_totalf) ' -- diffCN0_Bq = ' num2str(CN0_B_total-CN0_B_totalq) ' -- diffCN0_Bfq = ' num2str(CN0_B_total-CN0_B_totalfq)]);





