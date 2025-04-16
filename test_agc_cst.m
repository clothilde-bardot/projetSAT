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

M = 1/Ti; % nmb de sorties du corrélateurs sur 1s



% Paramètres AGC
v_ref = 3;          % Niveau de sortie désiré
alpha = 0.5;       % Coefficient du filtre de boucle
val_A = 1;           % Gain initial du VGA
val_Aq = val_A;




%BOUCLE DE SIMULATION

r_out = zeros(1,N);     % Signal de sortie

val_A_history = zeros(1,M); %historique gain AGC

r_out_history = zeros(1,M);

error_history = zeros(1,M);


for i = 1:nb_loops %nb d'itération
    disp(['Num Loops = ' num2str(i) '/' num2str(nb_loops)])

    for n = 1:M
        
        %AGC

        % Gain VGA
        r_in_gain = val_Aq*r_in;

        % Mesure de l'amplitude (detector)
        current_level = abs(r_in_gain(n));

        % Comparaison
        error = v_ref - current_level;

        % Signal de sortie
        r_out = r_in_gain;

        %MAJ du gain
        val_Aq = val_Aq + alpha * error;

        % Stockage pour visualisation
        val_A_history(n) = val_Aq;
        r_out_history(n) = mean(r_out);        
        error_history(n) = error;

    end 
end



%% Visualisation
figure;

% Signal d'entrée
subplot(3,1,1);
plot(t, r_in, 'b', 'LineWidth', 1.5);
title('Signal d''entrée x(t)');
xlabel('Temps (s)');
ylabel('Amplitude');
ylim([0 10]);
grid on;

% Signal de sortie
subplot(3,1,2);
plot(t, r_out, 'r', 'LineWidth', 1.5);
hold on;
yline(v_ref, '--k', 'v_{ref}');
title('Signal de sortie y(t) après AGC');
xlabel('Temps (s)');
ylabel('Amplitude');
ylim([0 10]);
grid on;

% Évolution du gain
subplot(3,1,3);
plot(t, val_A_history, 'g', 'LineWidth', 1.5);
title('Évolution du gain du VGA');
xlabel('Temps (s)');
ylabel('Gain');
grid on;
