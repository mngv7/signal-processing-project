%%%% EGB242 Assignment 2, Section 3 %%
% This file is a template for your MATLAB solution to Section 3.
%
% Before starting to write code, generate your data with the ??? as
% described in the assignment task.

%% Initialise workspace
clear all; close all; clc;
load DataA2 imagesReceived;

% Begin writing your MATLAB solution below this line.
%% Question 3.1
im1D = imagesReceived(1,:);

% Converting a received pixel stream to an image matrix
im2D = reshape(im1D, 480, 640) ;

% Displaying an image in a figure
figure;
imshow (im2D);

% Saving an image matrix as an image file 
imwrite (im2D , 'InitialImage.png') ;

%% Question 3.2
samples = length(imagesReceived);
f = 1000;
T = 1/1000;
% Creating the time vector
t = timevec(0, T * samples, samples); 
% Creating the Frequency vector
fs = freqvec(f, samples); 

im2D = im2D(1:samples);

figure;
plot(t, im2D);
xlim([0,307]);
ylim([-7,7]);
xlabel('Time [s]');
ylabel('Magnitdue');
title('Received image signal in Time Domain');

% IM2D = ft(im2D, fs);
% plot(fs, abs(IM2D));
% xlabel('Frequency [Hz]');
% ylabel('Magnitude');
% title('Received image signal in Frequency Domain');

%% Question 3.3

% Step input
step_input = ones(size(t));

% Define transfer function Filter 1
Pnum1 = [1]; 
Pden1 = [5.64e-5, 0.0167, 1]; 
Passive1tf = tf(Pnum1, Pden1);

% Define transfer function Filter 2
Pnum2 = [1.2e-3]; 
Pden2 = [5.64e-5, 5.9e-3, 1]; 
Passive2tf = tf(Pnum2, Pden2);

% Define transfer function Filter 1
Anum1 = [1.22e3]; 
Aden1 = [1, 2.44e3,1.44e6];
Active1tf = tf(Anum1, Aden1);

% Define transfer function Filter 2
Anum2 = [1]; 
Aden2 = [1, 2.44e3, 1.49e6]; 
Active2tf = tf(Anum2, Aden2);

%% LTI view
ltiview(Passive1tf);
ltiview(Passive2tf);
ltiview(Active1tf);
ltiview(Active2tf);

%%
% Use lsim to simulate the response of the feedback system
P1 = lsim(Passive1tf, step_input, t);
P2 = lsim(Passive2tf, step_input, t);
A1 = lsim(Active1tf, step_input, t);
A2 = lsim(Active2tf, step_input, t);

%% 3.3

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


