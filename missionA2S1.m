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
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Multiplexed Audio Signal (Time Domain)');

% Generate appropriate frequency vector.
samples = length(audioMultiplexNoisy);
f = freqvec(fs, samples);


% Convert audioMultiplexNoisy to frequency domain using Fourier transform.
AudioMultiplexNoisy =  ft(audioMultiplexNoisy, fs);

% Plot the frequency domain of audioMultiplexNoisy.
subplot(2, 1, 2);
plot(f, abs(AudioMultiplexNoisy), 'k');
title('Multiplexed Audio Signal (Frequency Domain)');
xlabel('Frequency [Hz]');
ylabel('Amplitude');

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
xlabel('Time (seconds)');
ylabel('Amplitude');
title('First Signal Demodulated (Time Domain)');

subplot(2, 1, 2);
plot(f, abs(FirstSignalDemodulated), 'k');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('First Signal Demodulated (Frequency Domain)');

%% Second signal centred at 24.33 kHz
carrierFrequency2 = 24330;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);

secondSignalDemodulated = audioMultiplexNoisy .* carrierSignal2;

%sound(secondSignalDemodulated, fs);

SecondSignalDemodulated =  ft(secondSignalDemodulated, fs);

figure;
subplot(2, 1, 1);
plot(t, secondSignalDemodulated, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Second Signal Demodulated (Time Domain)');

subplot(2, 1, 2);
plot(f, abs(SecondSignalDemodulated), 'k');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Second Signal Demodulated (Frequency Domain)');

%% Third signal centred at 40.27 kHz
carrierFrequency2 = 40270;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);

thirdSignalDemodulated = audioMultiplexNoisy .* carrierSignal2;

%sound(thirdSignalDemodulated, fs);

ThirdSignalDemodulated =  ft(thirdSignalDemodulated, fs);

figure;
subplot(2, 1, 1);
plot(t, thirdSignalDemodulated, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Third Signal Demodulated (Time Domain)');

subplot(2, 1, 2);
plot(f, abs(ThirdSignalDemodulated), 'k');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Third Signal Demodulated (Frequency Domain)');

%% Fourth signal centred at 56.17 kHz
carrierFrequency2 = 56170;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);

fourthSignalDemodulated = audioMultiplexNoisy .* carrierSignal2;

%sound(fourthSignalDemodulated, fs);

FourthSignalDemodulated =  ft(fourthSignalDemodulated, fs);

figure;
subplot(2, 1, 1);
plot(t, fourthSignalDemodulated, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Fourth Signal Demodulated (Time Domain)');

subplot(2, 1, 2);
plot(f, abs(FourthSignalDemodulated), 'k');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Fourth Signal Demodulated (Frequency Domain)');

%% Fifth signal centred at 72.29 kHz
carrierFrequency2 = 72290;

carrierSignal2 = cos(2*pi*carrierFrequency2*t);

fifthSignalDemodulated = audioMultiplexNoisy .* carrierSignal2;

%sound(fifthSignalDemodulated, fs);

FifthSignalDemodulated =  ft(fifthSignalDemodulated, fs);

figure;
subplot(2, 1, 1);
plot(t, fifthSignalDemodulated, 'k');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Fifth Signal Demodulated (Time Domain)');

subplot(2, 1, 2);
plot(f, abs(FifthSignalDemodulated), 'k');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Fifth Signal Demodulated (Frequency Domain)');
%% 1.3 Frequency and impulse response of the LTI system

% Vout(t)/Vin(t) = h(t)
% The output 'y' is the convolution of the input 'x' and the transfer function 'h'.

impulse = zeros(1, length(audioMultiplexNoisy));
impulse(1) = 1/(1/fs);
h_t = channel(sid, impulse, fs);
H_f = ft(h_t, fs);

figure;
plot(f, abs(H_f), 'r', f, abs(AudioMultiplexNoisy), 'b');
legend('H(f)', 'AudioMultiplexNoisy')
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Frequency response of H(f) against magntude spectrum of AudioMultiplexNoisy')
grid on;

%% 1.4 Reverse distortion

% Input(f) = Output(f) / H(f)
AudioMultiplexReverse = AudioMultiplexNoisy ./ H_f;
audioMultiplexReverse = ift(AudioMultiplexReverse, fs);

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

% meanValue1 = mean(demodAudioMultiplexReverse1);
% demodAudioMultiplexReverse1 = demodAudioMultiplexReverse1 - meanValue1;
figure;
subplot(2,1,1);
plot(t, demodAudioMultiplexReverse1, 'b')
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2);
plot(f, abs(DemodAudioMultiplexReverse1), 'r');
xlabel('Frequency [Hz]');
ylabel('Magnitude');


% fpass = 2e3;
% testfilter = lowpass(demodAudioMultiplexReverse1, fpass, fs);
% meanValue1 = mean(testfilter);
% testfilter = testfilter - meanValue1;
% figure;
% plot(t, testfilter);
% TestFilter = ft(testfilter, fs);
% figure;
% plot(f, abs(TestFilter));

% remove DC offset
meanValue1 = mean(demodAudioMultiplexReverse1);
demodAudioMultiplexReverse1 = demodAudioMultiplexReverse1 - meanValue1;

figure;
subplot(2,1,1);
plot(t, demodAudioMultiplexReverse1, 'b')
xlabel('Time (s)');
ylabel('Amplitude');
title('1st Audio time');
subplot(2,1,2);
plot(f, abs(DemodAudioMultiplexReverse1), 'r');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('1st Audio freq');

meanValue2 = mean(demodAudioMultiplexReverse2);
demodAudioMultiplexReverse2 = demodAudioMultiplexReverse2 - meanValue2;
figure;
subplot(2,1,1);
plot(t, demodAudioMultiplexReverse2, 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('2nd Audio time');
subplot(2,1,2);
plot(f, abs(DemodAudioMultiplexReverse2), 'r');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('2nd Audio freq');

meanValue3 = mean(demodAudioMultiplexReverse3);
demodAudioMultiplexReverse3 = demodAudioMultiplexReverse3 - meanValue3;
figure;
subplot(2,1,1);
plot(t, demodAudioMultiplexReverse3, 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('3rd Audio time');
subplot(2,1,2);
plot(f, abs(DemodAudioMultiplexReverse3), 'r');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('3rd Audio freq');

meanValue4 = mean(demodAudioMultiplexReverse4);
demodAudioMultiplexReverse4 = demodAudioMultiplexReverse4 - meanValue4;
figure;
subplot(2,1,1);
plot(t, demodAudioMultiplexReverse4, 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('4th Audio time');
subplot(2,1,2);
plot(f, abs(DemodAudioMultiplexReverse4), 'r');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('4th Audio freq');

meanValue5 = mean(demodAudioMultiplexReverse5);
demodAudioMultiplexReverse5 = demodAudioMultiplexReverse5 - meanValue5;
figure;
subplot(2,1,1);
plot(t, demodAudioMultiplexReverse5, 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('5th Audio time');
subplot(2,1,2);
plot(f, abs(DemodAudioMultiplexReverse5), 'r');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('5th Audio freq');

%% 1.5 Fully de-noising audio

signals = {demodAudioMultiplexReverse1, demodAudioMultiplexReverse2, demodAudioMultiplexReverse3, demodAudioMultiplexReverse4, demodAudioMultiplexReverse5};
newSignals = {'filteredAudioMultiplexReverse1', 'filteredAudioMultiplexReverse2', 'filteredAudioMultiplexReverse3', 'filteredAudioMultiplexReverse4', 'filteredAudioMultiplexReverse5'};
filteredSignals = struct();
fpass = 2e3;

for i = 1:length(signals)
    lowSignal = lowpass(signals{i}, fpass, fs, 'Steepness', 0.85, 'StopbandAttenuation', 60);
    filteredSignals.(newSignals{i}) = lowSignal;
end
filteredAudioMultiplexReverse1 = filteredSignals.filteredAudioMultiplexReverse1; filteredAudioMultiplexReverse2 = filteredSignals.filteredAudioMultiplexReverse2; filteredAudioMultiplexReverse3 = filteredSignals.filteredAudioMultiplexReverse3;
filteredAudioMultiplexReverse4 = filteredSignals.filteredAudioMultiplexReverse4; filteredAudioMultiplexReverse5 = filteredSignals.filteredAudioMultiplexReverse5;

%% 1.5 freq check plots
FilteredAudioMultiplexReverse1 = ft(filteredAudioMultiplexReverse1, fs); FilteredAudioMultiplexReverse2 = ft(filteredAudioMultiplexReverse2, fs); 
FilteredAudioMultiplexReverse3 = ft(filteredAudioMultiplexReverse3, fs); FilteredAudioMultiplexReverse4 = ft(filteredAudioMultiplexReverse4, fs); 
FilteredAudioMultiplexReverse5 = ft(filteredAudioMultiplexReverse5, fs); 

figure;
subplot(3,2,1);
plot(f, abs(FilteredAudioMultiplexReverse1));
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('First filtered audio in freq domain');

subplot(3,2,2);
plot(f, abs(FilteredAudioMultiplexReverse2));
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Second filtered audio in freq domain');

subplot(3,2,3);
plot(f, abs(FilteredAudioMultiplexReverse3));
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Third filtered audio in freq domain');

subplot(3,2,4);
plot(f, abs(FilteredAudioMultiplexReverse4));
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Fourth filtered audio in freq domain');

subplot(3,2,5);
plot(f, abs(FilteredAudioMultiplexReverse5));
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Fifth filtered audio in freq domain');





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
% Creates frequency vector suitable for plotting magnitude/phase spectrum
%
%  This function emulates the np.fft.freqvec function from python, but will
%  also make sure that the frequency vector has been shifted correctly, so
%  that the first index is for the lowest frequency, highest index is for
%  highest and that the middle frequency is DC.

%  https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html
%  
%  The frequency vector will be slightly different if the sequence is of 
%  even length or odd length. The specifics of this has to do with 
%  properties of the Discrete Fourier Transform (DTF), which is what the
%  FFT algorigthm actually computes. Will likely cover
%  this in your Digital signal processing class. The main thing to know is 
%  that this function will create a frequency vector that will ensure the DC
%  Component is at the correct location. For now, can just take our
%  word for it, and know what the function does and what it will return.
%
% 
%  For even length signals, our frequency vector will be of the form,
%        -fs/2 <= f < fs/2
%  For odd length signals, will be,
%        -fs/ 2 < f < fs/2
%
%  Args:
%  fs = sample frequency in Hz
%  n = length of the time vector/number of samples


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
