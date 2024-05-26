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
firstImage = imagesReceived(1,:);

numRows = 480;
numCols = 640;

% Converting a received pixel stream to an image matrix
reshapedFirstImage = reshape(firstImage, numRows, numCols);

% Displaying an image in a figure
figure;
imshow(reshapedFirstImage);

% Saving an image matrix as an image file 
imwrite(reshapedFirstImage, '1stImage.png');


%% Question 3.2

% Define the parameters
pixelRate = 1000; % pixels per second
numPixels = numel(firstImage); % returns the number of elements in an array (samples).
duration = numPixels / pixelRate; % total duration of the signal

% Construct the time vector
t = timevec(0, duration, numPixels);

% Construct the frequency vector
Fs = pixelRate; % Sampling frequency
f = freqvec(Fs, numPixels);

% Visualize the received signal in time domain
figure;
plot(t, firstImage);
xlim([0,307]);
ylim([-7,7]);
xlabel('Time [s]');
ylabel('Amplitude');
title('Received image signal in Time Domain');

% Visualize the received signal in frequency domain
FirstImage = ft(firstImage);
figure;
plot(f, abs(FirstImage));
title('Received Signal in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% Question 3.3

% Define transfer function Filter 1
Pnum1 = [1.2e-3]; 
Pden1 = [5.64e-5, 5.9e-3, 1]; 
Passive1tf = tf(Pnum1, Pden1);

% Define transfer function Filter 2
Pnum2 = [1]; 
Pden2 = [5.64e-5, 16.7e-3, 1]; 
Passive2tf = tf(Pnum2, Pden2);

% Define transfer function Filter 1
Anum1 = [1.49e6]; 
Aden1 = [1, 2.44e3,1.49e6];
Active1tf = tf(Anum1, Aden1);

% Define transfer function Filter 2
Anum2 = [1]; 
Aden2 = [1, 2.44e3, 1.49e6]; 
Active2tf = tf(Anum2, Aden2);

% %% LTI view
% ltiview(Passive1tf);
% ltiview(Passive2tf);
% ltiview(Active1tf);
% ltiview(Active2tf);

%% Use lsim to simulate the response of the feedback system
% P1 = lsim(Passive1tf, firstImage, t);
% P2 = lsim(Passive2tf, firstImage, t);
% A1 = lsim(Active1tf, firstImage, t);
% A2 = lsim(Active2tf, firstImage, t);

%% All images are best cleaned with Active filter 1.
%% 3.4
% First image

% Apply the filter to the received signal
 firstImage_filtered = lsim(Active1tf, firstImage, t);

% Reshape the filtered 1D signal back into a 2D image matrix
 reshapedFiltFirstImage = reshape(firstImage_filtered, numRows, numCols);

% Display the clean image
figure;
 imshow(reshapedFiltFirstImage);
 title('Filtered Landing Site Image');

% Visualize the clean signal in time domain
figure;
plot(t, firstImage_filtered);
xlim([0,307]);
ylim([-2,2]);
title('Filtered Signal in Time Domain');
xlabel('Time (s)');
ylabel('Pixel Intensity');

% Visualize the clean signal in frequency domain
FirstImage_filtered = ft(firstImage_filtered);
figure;
plot(f, abs(FirstImage_filtered));
title('Filtered Signal in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


%% 3.5
for i = 2:4
    image = imagesReceived(i, :);
    reshapedImage = reshape(image, numRows, numCols);
    imgString = [num2str(i) 'image.png'];
    imwrite(reshapedImage, imgString);

    figure;
    imshow(reshapedImage);

    image_filtered = lsim(Active1tf, image, t);

    reshapedFiltImage = reshape(image_filtered, numRows, numCols);

    figure;
    imshow(reshapedFiltImage);
    imgTitleString = ['Filtered Landing Site Image' num2str(i)];
    title(imgTitleString);

    figure;
    plot(t, image_filtered);
    xlim([0,307]);
    ylim([-2,2]);
    title('Filtered Signal in Time Domain');
    xlabel('Time (s)');
    ylabel('Pixel Intensity');

    Image_filteredP1 = ft(image_filtered);
    figure;
    plot(f, abs(Image_filteredP1));
    title('Filtered Signal in Frequency Domain');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
end
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

function fourierTransform = ft(freq)
% Applies complex Fourier transform with the assigned sampling frequency
%
%   Args:
%   freq = frequency/signal input
%   n = sample frequency in Hz

    fourierTransform = fftshift(fft(freq));
end
