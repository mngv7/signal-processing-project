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
ylabel('Amplitude');

% 1.2 Demodulate the signals in the frequency domain, listen, and plot.

% Signals probably need filtering?

% Frequencies to demodulate:
% 24.33 kHz
% 40.27 kHz
% 56.17 kHz
% 72.29 kHz

% First signal centred at 24.33 kHz.
carrierFrequency = 24330;

carrierSignal = cos(2*pi*carrierFrequency*t);

firstSignalDemodulated = audioMultiplexNoisy .* carrierSignal;

%sound(firstSignalDemodulated, fs);

FirstSignalDemodulated =  fftshift(fft(firstSignalDemodulated)) / fs;

figure;
subplot(2, 1, 1);
plot(t, firstSignalDemodulated, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('First Signal Demodulated (Time Domain)');

subplot(2, 1, 2);
plot(f, FirstSignalDemodulated, 'k');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('First Signal Demodulated (Frequency Domain)');

% Second signal centred at 40.27 kHz
carrierFrequency = 40270;

carrierSignal = cos(2*pi*carrierFrequency*t);

secondSignalDemodulated = audioMultiplexNoisy .* carrierSignal;

%sound(secondSignalDemodulated, fs);

SecondSignalDemodulated =  fftshift(fft(secondSignalDemodulated)) / fs;

figure;
subplot(2, 1, 1);
plot(t, secondSignalDemodulated, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Second Signal Demodulated (Time Domain)');

subplot(2, 1, 2);
plot(f, SecondSignalDemodulated, 'k');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Second Signal Demodulated (Frequency Domain)');

% Third signal centred at 56.17 kHz
carrierFrequency = 56170;

carrierSignal = cos(2*pi*carrierFrequency*t);

thirdSignalDemodulated = audioMultiplexNoisy .* carrierSignal;

%sound(thirdSignalDemodulated, fs);

ThirdSignalDemodulated =  fftshift(fft(thirdSignalDemodulated)) / fs;

figure;
subplot(2, 1, 1);
plot(t, thirdSignalDemodulated, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Third Signal Demodulated (Time Domain)');

subplot(2, 1, 2);
plot(f, ThirdSignalDemodulated, 'k');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Third Signal Demodulated (Frequency Domain)');

% Fourth signal centred at 72.29 kHz
carrierFrequency = 72290;

carrierSignal = cos(2*pi*carrierFrequency*t);

fourthSignalDemodulated = audioMultiplexNoisy .* carrierSignal;

%sound(fourthSignalDemodulated, fs);

FourthSignalDemodulated =  fftshift(fft(fourthSignalDemodulated)) / fs;

figure;
subplot(2, 1, 1);
plot(t, fourthSignalDemodulated, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Fourth Signal Demodulated (Time Domain)');

subplot(2, 1, 2);
plot(f, FourthSignalDemodulated, 'k');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Fourth Signal Demodulated (Frequency Domain)');

% 1.3 Frequency and impulse response of the LTI system

% Vout(t)/Vin(t) = h(t)
% The output 'y' is the convolution of the input 'x' and the transfer function 'h'.


y = channel(sid, x, fs);
Y =  fftshift(fft(y)) / fs;
