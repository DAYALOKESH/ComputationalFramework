% Okumura-Hata Model Path Loss Calculator
% This script implements the Okumura-Hata model for various environments.
% It calculates and plots Path Loss vs Distance.

clear;
clc;
close all;

fprintf('------------------------------------------------------------\n');
fprintf('Okumura-Hata Path Loss Model Simulation\n');
fprintf('------------------------------------------------------------\n');

% 1. User Input Section
% Defaults based on user request: f=970 MHz, h_t=52m, h_r=2.4m
default_f = 970;
default_ht = 52;
default_hr = 2.4;
default_dmax = 20;

fprintf('Enter parameters (press Enter to use defaults):\n');

% Frequency
f_input = input(sprintf('Frequency in MHz [Default: %.1f]: ', default_f));
if isempty(f_input)
    f = default_f;
else
    f = f_input;
end

% Transmitter Height
ht_input = input(sprintf('Transmitter Height (h_t) in meters [Default: %.1f]: ', default_ht));
if isempty(ht_input)
    h_t = default_ht;
else
    h_t = ht_input;
end

% Receiver Height
hr_input = input(sprintf('Receiver Height (h_r) in meters [Default: %.1f]: ', default_hr));
if isempty(hr_input)
    h_r = default_hr;
else
    h_r = hr_input;
end

% Max Distance
dmax_input = input(sprintf('Max Distance in km [Default: %.1f]: ', default_dmax));
if isempty(dmax_input)
    d_max = default_dmax;
else
    d_max = dmax_input;
end

% Basic Validation
if f < 150 || f > 1500
    warning('Frequency is outside the standard Hata model range (150-1500 MHz). Results may be inaccurate.');
end

% 2. Distance Vector
d = 0.1 : 0.1 : d_max; % Distance from 0.1 km to d_max km

% 3. Calculations

% Pre-calculate common terms
log_f = log10(f);
log_ht = log10(h_t);
% Common A, B, C, D terms for Urban
% L_urban = A + B*log_f - C*log_ht - a(h_r) + (D - E*log_ht)*log_d
A = 69.55;
B = 26.16;
C = 13.82;
D = 44.9;
E = 6.55;

% Correction factor a(h_r)

% Case 1: Urban (Small/Medium City)
% a(h_r) = (1.1*log(f) - 0.7)*h_r - (1.56*log(f) - 0.8)
a_hr_small_medium = (1.1 * log_f - 0.7) * h_r - (1.56 * log_f - 0.8);

L_urban_sm = A + B * log_f - C * log_ht - a_hr_small_medium + (D - E * log_ht) * log10(d);


% Case 2: Urban (Large City)
% a(h_r) depends on frequency
if f >= 300
    a_hr_large = 3.2 * (log10(11.75 * h_r))^2 - 4.97;
else
    a_hr_large = 8.29 * (log10(1.54 * h_r))^2 - 1.1;
end

L_urban_large = A + B * log_f - C * log_ht - a_hr_large + (D - E * log_ht) * log10(d);


% Case 3: Suburban
% Based on Urban (Small/Medium) usually
% L_suburban = L_urban - 2*(log(f/28))^2 - 5.4
L_suburban = L_urban_sm - 2 * (log10(f / 28))^2 - 5.4;


% Case 4: Rural (Open Area)
% Based on Urban (Small/Medium)
% L_rural = L_urban - 4.78*(log(f))^2 + 18.33*log(f) - 40.94
L_rural = L_urban_sm - 4.78 * (log_f)^2 + 18.33 * log_f - 40.94;


% 4. Plotting
figure;
plot(d, L_urban_large, 'r-', 'LineWidth', 2); hold on;
plot(d, L_urban_sm, 'b-', 'LineWidth', 2);
plot(d, L_suburban, 'g-', 'LineWidth', 2);
plot(d, L_rural, 'k-', 'LineWidth', 2);

grid on;
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
title(sprintf('Okumura-Hata Path Loss Model\n(f=%.1f MHz, Tx=%.1fm, Rx=%.1fm)', f, h_t, h_r));
legend('Urban (Large City)', 'Urban (Small/Medium City)', 'Suburban', 'Rural (Open Area)', 'Location', 'Best');

fprintf('\nCalculation Complete.\n');
fprintf('Displaying plot...\n');
