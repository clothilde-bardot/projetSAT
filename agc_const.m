%% Initialisation
clear all; close all; clc;

% Paramètres temporels
Fs = 1000;          % Fréquence d'échantillonnage (Hz)
T = 1;              % Durée totale (s)
t = 0:1/Fs:T-1/Fs;  % Vecteur temps
N = length(t);      % Nombre d'échantillons

% Paramètres AGC
Y_ref = 3;          % Niveau de sortie désiré (Y_{ref})
alpha = 0.01;       % Coefficient du filtre de boucle
gain = 1;           % Gain initial du VGA

% Signal d'entrée constant (x(t) = 5)
x = 5 * ones(1,N);  % Signal d'entrée

%% Boucle AGC
y = zeros(1,N);     % Signal de sortie
gain_history = zeros(1,N); % Historique du gain

for n = 1:N
    % Étape 1: VGA (amplification du signal)
    y(n) = gain * x(n);
    
    % Étape 2: Level detector (mesure de l'amplitude)
    current_level = abs(y(n)); % Détection d'amplitude simple
    
    % Étape 3: Calcul de l'erreur
    error = Y_ref - current_level;
    
    % Étape 4: Loop filter (ajustement progressif du gain)
    gain = gain + alpha * error;
    
    % Stockage pour visualisation
    gain_history(n) = gain;
end

%% Visualisation
figure;

% Signal d'entrée
subplot(3,1,1);
plot(t, x, 'b', 'LineWidth', 1.5);
title('Signal d''entrée x(t)');
xlabel('Temps (s)');
ylabel('Amplitude');
ylim([0 6]);
grid on;

% Signal de sortie
subplot(3,1,2);
plot(t, y, 'r', 'LineWidth', 1.5);
hold on;
yline(Y_ref, '--k', 'Y_{ref}');
title('Signal de sortie y(t) après AGC');
xlabel('Temps (s)');
ylabel('Amplitude');
ylim([0 6]);
grid on;

% Évolution du gain
subplot(3,1,3);
plot(t, gain_history, 'g', 'LineWidth', 1.5);
title('Évolution du gain du VGA');
xlabel('Temps (s)');
ylabel('Gain');
grid on;