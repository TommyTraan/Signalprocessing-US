%% 2D Linear imaging systems
% Exercise 7 in MEDT4165 
%
% Tommy Duy Tran

clear all; clc;


%% Part 1

showImages('Image.bmp', 1, 'Original line');
showImages('Image_translated.bmp', 2, 'Line translated');
showImages('Image_rotated1.bmp', 3, 'Line rotated 45 degrees');
showImages('Image_rotated2.bmp', 4, 'Line rotated 90 degrees');
showImages('Image_thick.bmp', 5, 'Line with increased line width');
showImages('Image_long.bmp', 6, 'Line with increased line length');


%% Part 2

im = showImages('lena.jpeg', 7, 'Lena');

showFiltering(im, 1, 1, 8)
showFiltering(im, 5, 5, 9)
showFiltering(im, 15, 15, 10)
showFiltering(im, 30, 30, 11)


%% Functions

% Part 1
function [PdB, xfreq, yfreq] = calculateFFT2(s, Nfft)
% Calculate 2D FFT of signal, returns signal of transform in decibel
    dim = 30E-2;        % cm
    dx = dim/Nfft;
    dy = dim/Nfft;
    
    xfreq = linspace(-0.5, 0.5, Nfft+1)*(1/dx); xfreq(end) = [];
    yfreq = linspace(-0.5, 0.5, Nfft+1)*(1/dy); yfreq(end) = [];
    
    PdB = 20*log10(abs(fftshift(fft2(s, Nfft, Nfft))));   % 20log10(A)
end


function im = showImages(fileName, figNumber, figTitle)
% Shows the image and 2D FFT of image
    
    % Open file
    im = imread(fileName);
    im = im(:,:,1);

    % Set figure size
    f = figure(figNumber);
    f.Position(3) = f.Position(3)*1.5;

    subplot(1,2,1)
    imagesc(im), axis('image'), xlabel('Length [mm]')
    title('Image in spatial domain')
    

    [PdB, xfreq, yfreq] = calculateFFT2(im, size(im,1));
    
    dyn = 50;        % dynamic range (contrast)
    gain = -100;     % gain (brightness)
    
    subplot(1,2,2)
    imagesc(xfreq, yfreq, PdB)
    axis('image'), clim([-dyn 0]-gain), xlabel('Frequency [Hz]')
    title('Image in Fourier domain')

    colormap(gray)
    sgtitle(figTitle)  % grid title
end


% Part 2
function PSF = lpFilter(Nx, Ny)
% Creates filter coeffictients of a lowpass filter
    PSF = ones(Nx,Ny)/(Nx*Ny);
end


function showFiltering(s, Nx, Ny, figNumber)
% Shows the filter, 2D FFT of filtered image and the filtered image
    % Create filter
    PSF = lpFilter(Nx, Ny);
    
    % Power spectra: 2D FFT of the filter matrix
    [PdB_filter, xfreq, yfreq] = calculateFFT2(PSF, size(s,1));

    % Resulting image
    im_filtered = filter2(PSF, s);

    % Power spectra: 2D FFT of filtered image
    [PdB_im_filtered, ~, ~] = calculateFFT2(im_filtered, size(s,1));
    
    % Create figure
    f = figure(figNumber);
    f.Position(3:4) = [f.Position(3)*2 f.Position(4)*0.8];

    subplot(1,3,1)
    imagesc(PdB_filter)
    axis('image'), xlabel('Frequency [Hz]')
    title('Filter matrix in Fourier domain')    
    
    subplot(1,3,2)
    imagesc(xfreq, yfreq, PdB_im_filtered)
    axis('image'), xlabel('Frequency [Hz]')
    title('Filtered image in Fourier domain')

    subplot(1,3,3)
    imagesc(im_filtered)
    axis('image'), xlabel('Length [mm]')
    title('Resulting filtered image in spatial domain')
    
    colormap(gray)
    sgtitle(sprintf('2D low pass filtering: Nx = %d, Ny = %d', Nx, Ny))
end