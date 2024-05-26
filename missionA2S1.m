%% EGB242 Assignment 2, Section 1 %%
% This file is a template for your MATLAB solution to Section 1.
%
% Before starting to write code, generate your data with the ??? as
% described in the assignment task.

%% Initialise workspace
clear all; %close all;
load DataA2 audioMultiplexNoisy fs sid;

% Begin writing your MATLAB solution below this line.

%% 1.1 Plot the time and frequency domain of audioMultiplexNoisy

% Generate appropriate time vector.
samples = length(audioMultiplexNoisy);
t = timevec(0, samples / fs, samples);


% Plot the time domain of audioMultiplexNoisy.
figure;
subplot(2, 1, 1);
plot(t, audioMultiplexNoisy, 'k');
xlabel('Time (seconds)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('Multiplexed Audio Signal (Time Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

% Generate appropriate frequency vector.
samples = length(audioMultiplexNoisy);
f = freqvec(fs, samples);


% Convert audioMultiplexNoisy to frequency domain using Fourier transform.
AudioMultiplexNoisy =  ft(audioMultiplexNoisy, fs);

% Plot the frequency domain of audioMultiplexNoisy.
subplot(2, 1, 2);
plot(f, abs(AudioMultiplexNoisy), 'k');
title('Multiplexed Audio Signal (Frequency Domain)', 'FontSize', 30);
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
%% 1.2 Demodulate the signals in the frequency domain, listen, and plot.

% Signals probably need filtering?

% Frequencies to demodulate:
% 8.05  kHz
% 24.33 kHz
% 40.27 kHz
% 56.17 kHz
% 72.29 kHz

% First signal centred at 8.05 kHz.
carrierFrequency2 = 8050;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);

firstSignalDemodulated = audioMultiplexNoisy .* carrierSignal2;

%sound(firstSignalDemodulated, fs);

FirstSignalDemodulated =  ft(firstSignalDemodulated, fs);

figure;
subplot(2, 1, 1);
plot(t, firstSignalDemodulated, 'k');
xlabel('Time (seconds)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('First Signal Demodulated (Time Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(2, 1, 2);
plot(f, abs(FirstSignalDemodulated), 'k');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('First Signal Demodulated (Frequency Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

%% Second signal centred at 24.33 kHz
carrierFrequency2 = 24330;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);

secondSignalDemodulated = audioMultiplexNoisy .* carrierSignal2;

%sound(secondSignalDemodulated, fs);

SecondSignalDemodulated =  ft(secondSignalDemodulated, fs);

figure;
subplot(2, 1, 1);
plot(t, secondSignalDemodulated, 'k');
xlabel('Time (seconds)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('Second Signal Demodulated (Time Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(2, 1, 2);
plot(f, abs(SecondSignalDemodulated), 'k');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('Second Signal Demodulated (Frequency Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

%% Third signal centred at 40.27 kHz
carrierFrequency2 = 40270;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);

thirdSignalDemodulated = audioMultiplexNoisy .* carrierSignal2;

%sound(thirdSignalDemodulated, fs);

ThirdSignalDemodulated =  ft(thirdSignalDemodulated, fs);

figure;
subplot(2, 1, 1);
plot(t, thirdSignalDemodulated, 'k');
xlabel('Time (seconds)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('Third Signal Demodulated (Time Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(2, 1, 2);
plot(f, abs(ThirdSignalDemodulated), 'k');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('Third Signal Demodulated (Frequency Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

%% Fourth signal centred at 56.17 kHz
carrierFrequency2 = 56170;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);

fourthSignalDemodulated = audioMultiplexNoisy .* carrierSignal2;

%sound(fourthSignalDemodulated, fs);

FourthSignalDemodulated =  ft(fourthSignalDemodulated, fs);

figure;
subplot(2, 1, 1);
plot(t, fourthSignalDemodulated, 'k');
xlabel('Time (seconds)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('Fourth Signal Demodulated (Time Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(2, 1, 2);
plot(f, abs(FourthSignalDemodulated), 'k');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('Fourth Signal Demodulated (Frequency Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

%% Fifth signal centred at 72.29 kHz
carrierFrequency2 = 72290;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);

fifthSignalDemodulated = audioMultiplexNoisy .* carrierSignal2;

%sound(fifthSignalDemodulated, fs);

FifthSignalDemodulated =  ft(fifthSignalDemodulated, fs);

figure;
subplot(2, 1, 1);
plot(t, fifthSignalDemodulated, 'k');
xlabel('Time (seconds)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('Fifth Signal Demodulated (Time Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(2, 1, 2);
plot(f, abs(FifthSignalDemodulated), 'k');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('Fifth Signal Demodulated (Frequency Domain)', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

%% 1.3 Frequency and impulse response of the LTI system

% Vout(t)/Vin(t) = h(t)
% The output 'y' is the convolution of the input 'x' and the transfer function 'h'.

impulse = zeros(1, length(audioMultiplexNoisy));
impulse(1) = 1/(1/fs);
h_t = channel(sid, impulse, fs);
H_f = ft(h_t, fs);
figure;
plot(f, abs(H_f), 'r', f, abs(AudioMultiplexNoisy), 'k');
legend('H(f)', 'AudioMultiplexNoisy')
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('Frequency response of H(f) against magntude spectrum of AudioMultiplexNoisy', 'FontSize', 30)
ax = gca;
ax.FontSize = 20;

grid on;

%% 1.4 Reverse distortion

% Input(f) = Output(f) / H(f)
AudioMultiplexReverse = AudioMultiplexNoisy ./ H_f;
audioMultiplexReverse = ift(AudioMultiplexReverse, fs);

figure;
plot(f, abs(AudioMultiplexReverse));
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('Frequency spectrum of Audio Multiplex Signal with reversed distortion', 'FontSize',30);
ax = gca;
ax.FontSize = 20;

% Demodulate signal again
carrierFrequency1 = 8050;

carrierSignal1 = cos(2*pi*carrierFrequency1*t);
demodAudioMultiplexReverse1 = audioMultiplexReverse .* carrierSignal1;

DemodAudioMultiplexReverse1 = ft(demodAudioMultiplexReverse1, fs);

%sound(demodAudioMultiplexReverse1, fs);

%%
% Demodulate signal again
carrierFrequency2 = 24330;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);
demodAudioMultiplexReverse2 = audioMultiplexReverse .* carrierSignal2;

DemodAudioMultiplexReverse2 = ft(demodAudioMultiplexReverse2, fs);

%sound(DemodAudioMultiplexReverse2, fs);

%%
% Demodulate signal again
carrierFrequency3 = 40270;

carrierSignal3 = cos(2*pi*carrierFrequency3*t);
demodAudioMultiplexReverse3 = audioMultiplexReverse .* carrierSignal3;

DemodAudioMultiplexReverse3 = ft(demodAudioMultiplexReverse3, fs);

%sound(DemodAudioMultiplexReverse3, fs);

%%
% Demodulate signal again
carrierFrequency4 = 56170;

carrierSignal4 = cos(2*pi*carrierFrequency4*t);
demodAudioMultiplexReverse4 = audioMultiplexReverse .* carrierSignal4;

DemodAudioMultiplexReverse4 = ft(demodAudioMultiplexReverse4, fs);

%sound(DemodAudioMultiplexReverse4, fs);

%%
% Demodulate signal again
carrierFrequency5 = 72290;

carrierSignal5 = cos(2*pi*carrierFrequency5*t);
demodAudioMultiplexReverse5 = audioMultiplexReverse .* carrierSignal5;

DemodAudioMultiplexReverse5 = ft(demodAudioMultiplexReverse5, fs);

%sound(DemodAudioMultiplexReverse5, fs);

%% 1.4 Plots

% plot in time and frequency domains


% remove DC offset
meanValue1 = mean(demodAudioMultiplexReverse1);
demodAudioMultiplexReverse1 = demodAudioMultiplexReverse1 - meanValue1;


subplot(5,2,1);
plot(t, demodAudioMultiplexReverse1, 'b')
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('1st Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,2);
plot(f, abs(DemodAudioMultiplexReverse1), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('1st Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

meanValue2 = mean(demodAudioMultiplexReverse2);
demodAudioMultiplexReverse2 = demodAudioMultiplexReverse2 - meanValue2;
subplot(5,2,3);
plot(t, demodAudioMultiplexReverse2, 'b');
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('2nd Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,4);
plot(f, abs(DemodAudioMultiplexReverse2), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('2nd Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

meanValue3 = mean(demodAudioMultiplexReverse3);
demodAudioMultiplexReverse3 = demodAudioMultiplexReverse3 - meanValue3;
subplot(5,2,5);
plot(t, demodAudioMultiplexReverse3, 'b');
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('3rd Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,6);
plot(f, abs(DemodAudioMultiplexReverse3), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('3rd Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

meanValue4 = mean(demodAudioMultiplexReverse4);
demodAudioMultiplexReverse4 = demodAudioMultiplexReverse4 - meanValue4;
subplot(5,2,7);
plot(t, demodAudioMultiplexReverse4, 'b');
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('4th Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,8);
plot(f, abs(DemodAudioMultiplexReverse4), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('4th Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

meanValue5 = mean(demodAudioMultiplexReverse5);
demodAudioMultiplexReverse5 = demodAudioMultiplexReverse5 - meanValue5;
subplot(5,2,9);
plot(t, demodAudioMultiplexReverse5, 'b');
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('5th Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,10);
plot(f, abs(DemodAudioMultiplexReverse5), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('5th Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

%% freq spectrums together
figure;
subplot(5,1,1);
plot(f, abs(DemodAudioMultiplexReverse1), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('1st Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(5,1,2);
plot(f, abs(DemodAudioMultiplexReverse2), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('2nd Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(5,1,3);
plot(f, abs(DemodAudioMultiplexReverse3), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('3rd Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(5,1,4);
plot(f, abs(DemodAudioMultiplexReverse4), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('4th Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(5,1,5);
plot(f, abs(DemodAudioMultiplexReverse5), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('5th Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

%% 1.5 Fully de-noising audio

signals = {demodAudioMultiplexReverse1, demodAudioMultiplexReverse2, demodAudioMultiplexReverse3, demodAudioMultiplexReverse4, demodAudioMultiplexReverse5};
newSignals = {'filteredAudioMultiplexReverse1', 'filteredAudioMultiplexReverse2', 'filteredAudioMultiplexReverse3', 'filteredAudioMultiplexReverse4', 'filteredAudioMultiplexReverse5'};
filteredSignals = struct();
fpass = 2e3;
bandfpass = [1.8e3 3e3];

for i = 1:length(signals)
    lowSignal = lowpass(signals{i}, fpass, fs, "Steepness", 0.9, "StopbandAttenuation", 60);
    bandSignal = bandstop(lowSignal, bandfpass, fs);
    filteredSignals.(newSignals{i}) = bandSignal;
end
filteredAudioMultiplexReverse1 = filteredSignals.filteredAudioMultiplexReverse1; filteredAudioMultiplexReverse2 = filteredSignals.filteredAudioMultiplexReverse2; filteredAudioMultiplexReverse3 = filteredSignals.filteredAudioMultiplexReverse3;
filteredAudioMultiplexReverse4 = filteredSignals.filteredAudioMultiplexReverse4; filteredAudioMultiplexReverse5 = filteredSignals.filteredAudioMultiplexReverse5;
%%
d = designfilt('lowpassiir', 'PassbandFrequency', fpass/(fs/2), ...
    'StopbandFrequency', (fpass+fpass*steepness)/(fs/2), ...
    'StopbandAttenuation', stopband, 'DesignMethod', 'butter');

% Visualize the filter using fvtool
h = fvtool(d);

% Set the analysis to magnitude response
h.Analysis = 'magnitude';

% Obtain the handle to the axes
ax = get(h, 'CurrentAxes');

% Set the Y-axis limits to focus on the range from 0 dB to -80 dB
ax.YLim = [-80 0];

% Optionally, limit the X-axis range, for example, from 0 to 0.5 (normalized frequency)
ax.XLim = [0 0.5];
%%
d = designfilt('bandstopiir', 'FilterOrder', 10, ...
    'HalfPowerFrequency1', bandfpass(1), 'HalfPowerFrequency2', bandfpass(2), ...
    'DesignMethod', 'butter', 'SampleRate', fs);

% Visualize the bandstop filter using fvtool
h = fvtool(d);

% Set the analysis to magnitude response
h.Analysis = 'magnitude';

% Obtain the handle to the axes
ax = get(h, 'CurrentAxes');

% Set the Y-axis limits to focus on the range from 0 dB to -80 dB
ax.YLim = [-80 0];

% Optionally, limit the X-axis range to a specific range (e.g., from 0 to 0.5 normalized frequency)
ax.XLim = [0 0.5];
%%
% Design the lowpass filter
lowpassFilter = designfilt('lowpassiir', 'PassbandFrequency', fpass/(fs/2), ...
    'StopbandFrequency', (fpass + fpass * steepness)/(fs/2), ...
    'StopbandAttenuation', stopband, 'DesignMethod', 'butter');

% Design the bandstop filter
bandstopFilter = designfilt('bandstopiir', 'FilterOrder', 10, ...
    'HalfPowerFrequency1', bandfpass(1)/(fs/2), 'HalfPowerFrequency2', bandfpass(2)/(fs/2), ...
    'DesignMethod', 'butter');

% Convert the digitalFilter objects to second-order sections (SOS)
lowpassSOS = ss2sos(lowpassFilter);
bandstopSOS = ss2sos(bandstopFilter);

% Create dfilt objects from SOS
lowpassDfilt = dfilt.df2sos(lowpassSOS);
bandstopDfilt = dfilt.df2sos(bandstopSOS);

% Combine the filters
combinedFilter = dfilt.cascade(lowpassDfilt, bandstopDfilt);

% Visualize the combined filter using fvtool
h = fvtool(combinedFilter);

% Set the analysis to magnitude response
h.Analysis = 'magnitude';

% Obtain the handle to the axes
ax = get(h, 'CurrentAxes');

% Set the Y-axis limits to focus on the range from 0 dB to -80 dB
ax.YLim = [-80 0];

% Optionally, limit the X-axis range to a specific range (e.g., from 0 to 0.5 normalized frequency)
ax.XLim = [0 0.5];

% Function to convert digitalFilter object to SOS matrix
function sos = ss2sos(d)
    sos = d.Coefficients;
end
%% 1.5 freq check plots
FilteredAudioMultiplexReverse1 = ft(filteredAudioMultiplexReverse1, fs); FilteredAudioMultiplexReverse2 = ft(filteredAudioMultiplexReverse2, fs); 
FilteredAudioMultiplexReverse3 = ft(filteredAudioMultiplexReverse3, fs); FilteredAudioMultiplexReverse4 = ft(filteredAudioMultiplexReverse4, fs); 
FilteredAudioMultiplexReverse5 = ft(filteredAudioMultiplexReverse5, fs); 

subplot(5,2,1);
plot(t, filteredAudioMultiplexReverse1, 'b')
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('1st Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,2);
plot(f, abs(FilteredAudioMultiplexReverse1), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('1st Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(5,2,3);
plot(t, filteredAudioMultiplexReverse2, 'b');
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('2nd Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,4);
plot(f, abs(FilteredAudioMultiplexReverse2), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('2nd Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(5,2,5);
plot(t, filteredAudioMultiplexReverse3, 'b');
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('3rd Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,6);
plot(f, abs(FilteredAudioMultiplexReverse3), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('3rd Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(5,2,7);
plot(t, filteredAudioMultiplexReverse4, 'b');
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('4th Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,8);
plot(f, abs(FilteredAudioMultiplexReverse4), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('4th Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;

subplot(5,2,9);
plot(t, filteredAudioMultiplexReverse5, 'b');
xlabel('Time (s)', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
title('5th Audio time', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;
subplot(5,2,10);
plot(f, abs(FilteredAudioMultiplexReverse5), 'r');
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Magnitude', 'FontSize', 30);
title('5th Audio freq', 'FontSize', 30);
ax = gca;
ax.FontSize = 20;





%% helper functions
% function definitions in matlab either need to be in their own file,
% or can be in at the bottom of a script.
% if you want to these functions outside this lab, feel free to 
% move them into their own file, just make sure the filename is the same 
% as the function name, ie timevec.m, freqvec.m


function t=timevec(t0, t0_plus_T, n)
% Creates time vector, where upper limit is non-inclusive
%          t0 <= t < t0_plus_T
%   It is the responsibility of the user to ensure that for the use-case
%   that they want the lower limit included, and the upper-limit 
%   not included.
%
%   Args:
%   t0 = start time
%   t0_plus_T = end time (t0 + T)
%   n = number of samples

    t = linspace(t0, t0_plus_T, n + 1);
    t = t(1:end - 1);
end


function fourierTransform = ft(freq, n)
% Applies complex Fourier transform with the assigned sampling frequency
%
%   Args:
%   freq = frequency/signal input
%   n = sample frequency in Hz

    fourierTransform = fftshift(fft(freq)) / n;
end

function inversefourierTransform = ift(freq, n)
% Applies inverse Fourier transform with the assigned sampling frequency
%
%   Args:
%   freq = frequency/signal input
%   n = sample frequency in Hz

    inversefourierTransform = ifft(ifftshift(freq)) * n;
end


function f=freqvec(fs, n)



    % if is an even sequence length, generating the frequency vector
    % is just like doing it for our time vector in the timevec function
    if mod(n, 2) == 0
        f_str = sprintf('Generating freq vec\n [%.2f, %.2f)\n', -fs/2, fs/2);
        disp(f_str);
        % compute the frequency vector
        f = linspace(-fs / 2, fs / 2, n + 1);
        f = f(1:end - 1); 
    % otherwise is of odd length
    else        
        f_str = sprintf('Generating freq vec\n (%.2f, %.2f)\n', ...
            -(n -1)/2 * fs / n, (n -1)/2 * fs / n);
        disp(f_str);
        % compute the frequency vector
        f = linspace(-(n -1)/ 2, (n -1) / 2, n) * fs / n;
    end
end
