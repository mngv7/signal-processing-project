%% EGB242 Assignment 2, Section 1 %%
% This file is a template for your MATLAB solution to Section 1.
%
% Before starting to write code, generate your data with the ??? as
% described in the assignment task.

%% Initialise workspace
clear all; close all;
load DataA2 audioMultiplexNoisy fs sid;

% Begin writing your MATLAB solution below this line.

% 1.1 Plot the time and frequency domain of audioMultiplexNoisy

% Generate appropriate time vector.
samples = length(audioMultiplexNoisy);
t = linspace(0, samples / fs, samples + 1);
t(end) = [];

% Plot the time domain of audioMultiplexNoisy.
figure;
subplot(2, 1, 1);
plot(t, audioMultiplexNoisy, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Multiplexed Audio Signal (Time Domain)');

% Generate appropriate frequency vector.
samples = length(audioMultiplexNoisy);
f = linspace(-fs/2, fs/2, samples+1);
f(end) = [];

% Convert audioMultiplexNoisy to frequency domain using Fourier transform.
AudioMultiplexNoisy =  fftshift(fft(audioMultiplexNoisy)) / fs;

% Plot the frequency domain of audioMultiplexNoisy.
subplot(2, 1, 2);
plot(f, AudioMultiplexNoisy, 'k');
title('Multiplexed Audio Signal (Frequency Domain)');
xlabel('Frequency [Hz]');
ylabel('Magnitude');

% 1.2 Demodulate the signals in the frequency domain, listen, and plot.

% Signals probably need filtering?
% Plot time and frequency domain of each signal.

% Frequencies to demodulate:
% 24.33 kHz
% 40.27 kHz
% 56.17 kHz
% 72.29 kHz

% First signal centred at 24.33 kHz.
carrierFrequency = 24330;

carrierSignal = cos(2*pi*carrierFrequency*t);

firstSignalDemodulated = audioMultiplexNoisy .* carrierSignal;

%sound(firstSignalFiltered, fs);

% Second signal centred at 40.27 kHz
carrierFrequency = 40270;

carrierSignal = cos(2*pi*carrierFrequency*t);

secondSignalDemodulated = audioMultiplexNoisy .* carrierSignal;

%sound(secondSignalDemodulated, fs);

% Third signal centred at 56.17 kHz
carrierFrequency = 56170;

carrierSignal = cos(2*pi*carrierFrequency*t);

thirdSignalDemodulated = audioMultiplexNoisy .* carrierSignal;

%sound(thirdSignalDemodulated, fs);

% Fourth signal centred at 72.29 kHz
carrierFrequency = 72290;

carrierSignal = cos(2*pi*carrierFrequency*t);

fourthSignalDemodulated = audioMultiplexNoisy .* carrierSignal;

%sound(fourthSignalDemodulated, fs);

