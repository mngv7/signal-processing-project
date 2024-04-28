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

samples = length(audioMultiplexNoisy);
t = linspace(0, samples / fs, samples + 1);
t(end) = [];

% Plot original audio signal
figure;
subplot(2, 1, 1);
plot(t, audioMultiplexNoisy, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Multiplexed Audio Signal (Time Domain)');

samples = length(audioMultiplexNoisy);
f = linspace(-fs/2, fs/2, samples+1);
f(end) = [];

AudioMultiplexNoisy =  fftshift(fft(audioMultiplexNoisy)) / fs;

% Plot magnitude spectrum
subplot(2, 1, 2);
plot(f, AudioMultiplexNoisy, 'k');
title('Multiplexed Audio Signal (Frequency Domain)');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
