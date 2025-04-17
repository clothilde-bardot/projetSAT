%INITIALISATION
clear all;
close all;
clc;


Fs = 1000;   %fréquence d'échantillonage
T = 1;              % Durée totale (s)
t = 0:1/Fs:T-1/Fs;  % Vecteur temps
N = length(t);      % Nombre d'échantillons


r_in = 5* ones(1,N);     %entrée constante


nb_loops = 10;

Ti = 1e-3;

M = 1/Ti; % nmb de sorties du corrélateur sur 1s

total_iter = nb_loops * M;

% Paramètres AGC
v_ref = 3;          % Niveau de sortie désiré
K_0 = 0.1;       % Coefficient du filtre de boucle
g0 = 1;           % Gain initial du VGA
gq = g0;




%BOUCLE DE SIMULATION

r_out = zeros(1,N);     % Signal de sortie

g_history = zeros(1, total_iter);
r_out_history = zeros(1, total_iter);
error_history = zeros(1, total_iter);

k=1;

for i = 1:nb_loops %nb d'itération
    disp(['Num Loops = ' num2str(i) '/' num2str(nb_loops)])

    for n = 1:M
        %AGC

        % Gain VGA
        r_in_gain = gq*r_in;

        % Mesure de l'amplitude (detector)
        current_level = abs(r_in_gain(n));

        % Comparaison
        error = v_ref - current_level;

        % Signal de sortie
        r_out = r_in_gain(n);
        disp(['Itération ' num2str(n) ' : r\_out = ' num2str(r_out)]);
        %MAJ du gain
        %gq = gq + K_0 * error;
        gq = min(max(gq + K_0 * error, 0), 10);

        % Stockage pour visualisation
        g_history(k) = gq;
        r_out_history(k) = r_out;        
        error_history(k) = error;
        k=k+1;
    end 
end

t_agc = (0:total_iter-1) * Ti;

%% Visualisation
figure;

% Signal d'entrée
subplot(3,1,1);
plot(t, r_in, 'b', 'LineWidth', 1.5);
title('Signal d''entrée');
xlabel('Temps (s)');
ylabel('Amplitude');
ylim([0 10]);
grid on;

% Signal de sortie
subplot(3,1,2);
plot(t_agc, r_out_history, 'r', 'LineWidth', 1.5);
hold on;
yline(v_ref, '--k', 'v_{ref}');
title('Signal de sortie après AGC');
xlabel('Temps (s)');
ylabel('Amplitude');
ylim([0 10]);
grid on;

% Évolution du gain
subplot(3,1,3);
plot(t_agc, g_history, 'g', 'LineWidth', 1.5);
title('Évolution du gain');
xlabel('Temps (s)');
ylabel('Gain');
grid on;
