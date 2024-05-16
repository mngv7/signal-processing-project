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
imshow (im2D) ;

% Saving an image matrix as an image file 
imwrite (im2D , 'InitialImage.png') ;

%% Question 3.2
samples = length(im1D);
f = 1000;
T = 1/1000;
% Creating the time vector
t = timevec(0, T * samples, samples); 
% Creating the Frequency vector
fs = freqvec(f, samples); 

figure;
plot(t, im1D);
xlim([0,307]);
ylim([-7,7]);
xlabel('Time [s]');
ylabel('Magnitdue');
title('Received image signal in Time Domain');
%%
% IM1D = ft(im1D, fs);
% plot(fs, abs(IM1D));
% xlabel('Frequency [Hz]');
% ylabel('Magnitude');
% title('Received image signal in Frequency Domain');
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


