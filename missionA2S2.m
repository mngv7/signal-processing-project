%% EGB242 Assignment 2, Section 2 %%
% This file is a template for your MATLAB solution to Section 2.
%
% Before starting to write code, generate your data with the ??? as
% described in the assignment task.

%% Initialise workspace
clear all; close all;

%% 2.1

t2 = timevec(0, 20, 100000); % Time vector

% Step input
step_input = ones(size(t2));

% Step response
step_response = 2*t2 + 4*exp(-0.5*t2) - 4;

% Plotting
figure;
subplot(2,1,1);
plot(t2, step_input, 'b', 'LineWidth', 1.5); % Step input
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Step Input');
grid on;

subplot(2,1,2);
plot(t2, step_response, 'r', 'LineWidth', 1.5); % Step response
xlabel('Time (s)');
ylabel('Radians (rad)');
title('Step Response');
grid on;

%% 2.2

% Define transfer function G1
numerator_G1 = [1];
denominator_G1 = [1, 0.5, 1];
G1 = tf(numerator_G1, denominator_G1);

% Use lsim to simulate the response of the feedback system
y1 = lsim(G1, step_input, t2);

% Plot the step response
figure;
plot(t2, y1);
title('Step Response of Feedback System');
xlabel('Time (s)');
ylabel('Output');
grid on;

% This isn't good enough because it doesn't rotate all the way.

%% 2.3

% Underdamped
% Natural frequency: 1
% Damping ratio: 0.25

% Peak time: 3.24s
% Settling time: 16s
% Overshoot %: 44.4%

%% 2.4
    
% Define the values for Ggs and Hgs
Ggs = [0.1, 0.2, 0.5, 1, 2];
Hgs = [0.1, 0.2, 0.5, 1, 2];

% Loop over the values in Ggs while keeping Hgs constant
for i = 1:5
    Ggs_i = 1;
    num = Ggs_i;
    den = [1, 0.5, Hgs(i)];
    G2 = tf(num, den);
    
    y = lsim(G2, step_input, t2);
    
    subplot(2, 5, i); % Create a subplot grid (2 rows, 5 columns)
    plot(t2, y);
    title(['Step Response Kfb = ', num2str(Hgs(i))]);
    xlabel('Time (s)');
    ylabel('Output');
end

% Loop over the values in Hgs while keeping Ggs constant
for i = 1:5
    Hgs_i = 1;
    num2 = Hgs(i);
    den2 = [1, 0.5, Ggs(i)];
    G2 = tf(num2, den2);
    
    y = lsim(G2, step_input, t2);

    subplot(2, 5, 5+i); % Create a subplot grid (2 rows, 5 columns)
    plot(t2, y);
    title(['Step Response Kfwd = ', num2str(Ggs(i))]);
    xlabel('Time (s)');
    ylabel('Output');
end
sgtitle('System Responses for Different Ggs and Hgs Values');

%% 2.5

% Get the system parameters (2.3) of all the plots created in 2.4, using
% this determine which gain values (Hg(s) and Gg(s)) is suitable. Requires
% of the system:

% - Output must accurately travel between 0 to 2pi rads (so preferrably no
% oscillation).
% - The camera shouldn't rotate to quickly (peak time shouldn't be too
% fast, aim for 13 seconds).

% Define transfer function G1
Fnum = [1];
Fden = [1, 0.5, 0.158];
cameraTF = tf(Fnum, Fden);

% Use lsim to simulate the response of the feedback system
y2 = lsim(cameraTF, step_input, t2);

% Plot the step response
figure;
plot(t2, y2);
title('Step Response of Feedback System');
xlabel('Time (s)');
ylabel('Output');
grid on;

%% 2.6
startVoltage = 30/Fden; %Not sure yet
endVoltage = 210/Fden; %Not sure yet either
[startIm, finalIm] = cameraPan(startVoltage, endVoltage, cameraTF);
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
