%% Beamforming
% Exercise 6 in MEDT4165 
%
% Tommy Duy Tran

clear all; clc;


%% Load file and create variables
load simdata.mat;
c = 1540;                                   % m/s
channels = 1:1:size(RFdata,2);


%% Task 1: Time of Flight calculation

TOF = calculateTOF(c, 0.01, 0, elPosX);     % Calculate TOF
im = 20*log10(abs(RFdata));                 % Visualize data in log scale

figure(1)
imagesc(channels, RF_t, im), colormap(gray)
hold on
plot(channels, TOF,'linewidth', 1.5, 'Color', 'r'), legend('TOF')
xlabel('Channel number')
title('RF data on logarithmic scale, and time of flight')


%% Calculate TOF matrix

x = linspace(-2e-2, 2e-2, 256);            % x-coordinates
z = linspace(0, 4.5e-2, 512);              % z-coordinates
[X, Z, ElPosX] = meshgrid(x, z, elPosX);   % make X, Z and transducer element position matrices

% Matrix with TOF for each z-pos, x-pos and transducer element
TOF = calculateTOF(c, Z, X, ElPosX);


%% Task 2: Apodization

delayedData = interpTOF(RFdata, RF_t, TOF);     % interpolate channel data

bf_RFdata = squeeze(sum(delayedData,1));        % sum channel data and remove
                                                % the dimension of length 1

envelope = abs(hilbert(bf_RFdata));             % create envelope
envelope_dB = 20*log10(envelope);               % signal in dB

figure(2)
beamformedImage(channels, z, envelope_dB, 'Beamformed image (without apodization)')


%% Apodization: Hamming window

l = length(channels);
envelope_dB = apodization(delayedData, hamming(l));     % Apodization

figure(3)
beamformedImage(channels, z, envelope_dB, 'Beamformed image (apodization: Hamming window)')


%% Apodization: Hamming window on middle 32 elements, 0 on rest

l = 32;

% Creating apodization function
hamming_mid32 = zeros(length(channels),1);
hamming_mid32(length(hamming_mid32)/2 - l/2 + 1: length(hamming_mid32)/2 + l/2) = hamming(l);

envelope_dB = apodization(delayedData, hamming_mid32);      % Apodization

figure(4)
beamformedImage(channels, z, envelope_dB, 'Beamformed image (apodization: Hamming window on middle 32 elements)')


%% Apodization: Hamming window on first 32 elements, 0 on rest

% Creating apodization function
hamming_first32 = zeros(length(channels),1);
hamming_first32(1:l) = hamming(l);

envelope_dB = apodization(delayedData, hamming_first32);    % Apodization

figure(5)
beamformedImage(channels, z, envelope_dB, 'Beamformed image (apodization: Hamming window on first 32 elements)')


%% Task 3: Expanding aperture
xpos = 0;
fnumber = 1;

apod = generateApod(elPosX, xpos, z, fnumber);

figure(6)
imagesc(channels, z*1E3, transpose(apod))
xlabel('Channel number'), ylabel('Depth [mm]')
title('Apodization function, F# number = 1')


%% Imaging with expanding aperture

% Expanding aperture 
envelope_dB = expandingAperture(z, x, elPosX, delayedData, 1);
figure(7)
beamformedImage(channels, z, envelope_dB, 'Beamformed image (expanding aperture, F# = 1)')


%% Compare different F-numbers

figure(8)
subplot(1,2,1)
envelope_dB = expandingAperture(z, x, elPosX, delayedData, 0.1);
beamformedImage(channels, z, envelope_dB, 'F# = 0.1')

subplot(1,2,2)
envelope_dB = expandingAperture(z, x, elPosX, delayedData, 10);
beamformedImage(channels, z, envelope_dB, 'F# = 10')


%% In vivo data

load invivoData.mat;

% Calculate TOF
[X, Z, ElPosX] = meshgrid( x, z, elPosX); 
TOF = calculateTOF(c, Z, X, ElPosX);

% Interpolate data
delayedData = interpTOF(RFdata, RF_t, TOF); 

% Expanding aperture 
envelope_dB = expandingAperture(z, x, elPosX, delayedData, 1);

figure(9)
beamformedImage(channels, z, envelope_dB, 'Beamformed image, data in vivo, F# = 1')


%% Functions

function TOF = calculateTOF(c, z, x, x_e)
% Calculates time of flight (TOF)
% c: speed of sound
% (x,z): x and z [m] position of scatterer
% x_e: lateral position of elements in transducer
    t_tx = z/c;                          % time: pulse excitation to scatterer
    t_rx = sqrt((x-x_e).^2 + z.^2)./c;   % time: backscattered signal hits element
    TOF = t_tx + t_rx;
end


function beamformedImage(channels, z, envelope_dB, figTitle)
% Shows beamformed image of data
    dyn = 50;       % dynamic range (contrast)
    gain = -100;    % gain (brightness)

    imagesc(channels, z*1e3, envelope_dB)
    colormap('gray');
    clim([-dyn 0]-gain)
    xlabel('Channel number'), ylabel('Depth [mm]')
    title(figTitle)
end


function envelope_dB = apodization(delayedData, apod_func)
% Apodization on data
% Returns the envelope in dB for plotting
    apod = squeeze(sum(delayedData .* apod_func,1));
    envelope_dB = 20*log10(abs(hilbert(apod)));
end


function envelope_dB = expandingAperture(z, x, elPosX, delayedData, fnumber)
% Data generated with expanding aperture
% Returns envelope in dB for plotting
    apod_matrix = zeros(length(elPosX), length(z), length(x));
    for i = 1:length(x)
        apod_matrix(:,:,i) = generateApod(elPosX, x(i), z, fnumber);
    end

    envelope_dB = apodization(delayedData, apod_matrix);
end