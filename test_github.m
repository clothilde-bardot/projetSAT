clc; clear; close all;

% Paramètres du signal d'entrée
amplitude = 200e-6;  % 200 µV
f = 2e6;  % 2 MHz

% Fréquence d'échantillonnage
Fs = 16e6;

% Gain initial et bornes
max_gain = 6300;
min_gain = 1;  % Gain minimal pour éviter 0
gain = max_gain;  

% Niveau de sortie souhaité (Y_ref)
Y_ref = 50e-3;  % 50 mV

% Nombre d'échantillons et vecteur temps
N = 1e2;
dT = 1/Fs;
t = 0:dT:(N-1)*dT;

% Génération du bruit blanc
noise_power = (1e-9)^2 * (Fs/2);
noise = sqrt(noise_power) * randn(1,N);

% Signal d'entrée avec bruit
vin = amplitude * sin(2*pi*f*t) + noise;

% Initialisation des signaux
vin_amplified = zeros(1,N);
gain_values = zeros(1,N);
error_signal = zeros(1,N);

% Filtre de boucle : coefficient d'ajustement du gain
alpha = 0.1;

% Boucle AGC
for ii = 1:N
    % Amplification du signal
    vin_amplified(ii) = gain * vin(ii);
    
    % Calcul de l'erreur entre la sortie actuelle et la référence
    error_signal(ii) = Y_ref - abs(vin_amplified(ii));
    
    % Mise à jour du gain avec un filtre proportionnel
    gain = gain * (1 + alpha * error_signal(ii));
    
    % Saturation du gain pour éviter des valeurs extrêmes
    gain = max(min_gain, min(gain, max_gain));
    
    % Stockage du gain pour affichage
    gain_values(ii) = gain;
end

% Affichage des résultats
figure;
plot(t, vin_amplified);
title('Output Voltage after AGC');
xlabel('Time [s]');
ylabel('Amplitude [V]');
grid on;

figure;
plot(t, error_signal);
title('Error Signal');
xlabel('Time [s]');
ylabel('Error [V]');
grid on;

figure;
plot(t, gain_values);
title('Gain Evolution');
xlabel('Time [s]');
ylabel('Gain [V/V]');
grid on;
