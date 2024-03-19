%% Signal Processing Chain
% Exercise 5 in MEDT4165 
%
% Tommy Duy Tran

clear all; clc;



%% Task 1: Complex demodulation
[p, RF, f0, fs] = loadFile('phantomRFdata1.mat');

IQ = performTask1(RF, f0, fs, 1, 'Phantom Data 1: Power spectra before and after demodulation');



%% Task 2: Envelope detection
envelope = performTask2(p, RF, fs, IQ, 2, 'Phantom Data 1: Envelope and middle image line');



%% Task 3: Scanconversion and log-compression
performTask3(p, RF, envelope, 3, 'Phantom Data 1: Image');



%% Task 4: Settings
[p, RF, f0, fs] = loadFile('phantomRFdata2.mat');

IQ = performTask1(RF, f0, fs, 4, 'Phantom Data 2: Power spectra before and after demodulation');
envelope = performTask2(p, RF, fs, IQ, 5, 'Phantom Data 2: Envelope and middle image line');
performTask3(p, RF, envelope, 6, 'Phantom Data 2: Image');



%% Task 5: Fundamental vs. harmonic imaging
[p, RF, f0, fs] = loadFile('cardiacRFdata.mat');

IQ = performTask1(RF, f0, fs, 7, 'Cardiac Data: Power spectra before and after demodulation');
envelope = performTask2(p, RF, fs, IQ, 8, 'Cardiac Data: Envelope and middle image line');
performTask3(p, RF, envelope, 9, 'Cardiac Data: Image');



%% Functions
function [p, RF, f0, fs] = loadFile(fileName)
% Load the file
    load(fileName);
    RF = rf;
    f0 = p.f0_Hz;% * 1E-6;     % [MHz], transmit demodulation frequency
    fs = p.frs_Hz;% * 1E-6;    % [MHz], sample frequency
end


% Task 1a)
function dm = downmix(s, f0, fs)
% Downmixing the signal
    ts = 1/fs;                 % time step
    t = transpose((0:1:size(s,1)-1) * ts);   % time starting from 0 with timestep ts

    w = 2*pi*f0;               % omega
    dm = s .* exp(-1i*w*t);    % downmixing the signal
end

function lp = lpFilter(s, f0, fs)
% Lowpass filters the signal
    coeff = transpose(fir1(100, f0/(fs/2), 'low'));   % coefficients for filters
    lp = filter2(coeff, s);    % filtering the signal
end


% Task 1b)
function IQ = demodulation(s, f0, fs)
% Demodulation by downmixing and lowpass filtering
    dm = downmix(s, f0, fs);
    IQ = lpFilter(dm, f0, fs);
end


% Task 1c)
function PdB = calculatePower(s)
% Calculates the power spectra of signal
    Nfft = size(s,1);
    PdB = 20*log10(abs(fftshift(fft(s, Nfft))));   % 20log10(A)
end

function ave = averagePower(s)
% Averages the power over all image lines (rows)
    PdB = calculatePower(s);
    ave = mean(PdB,2);         % returns column vector containing mean of each row
end


% Task 2a)
function [depth, line] = getLine(p, RF, fs, index)
% Get wanted image line
    startdepth = p.startdepth_m;   % get data from struct
    c = 1540;            % [m/s]   % given parameter
    
    line = RF(:,index);
    
    ts = 1/fs;
    t = transpose((0:1:size(RF,1)-1) * ts);

    depth = t*c + startdepth;      % times * velocity = distance/depth, shift array with startdepth
end


% Task 2b)
function envelope = getEnvelope(IQ)
% Find envelope of signal
    envelope = 2 * abs(IQ);
end


% Task 3a)
function im = beamToCartesian(p, RF, envelope)
% Converts image from beam space to cartesian space
    theta0 = p.startangle_rad;          % get data from struct
    theta_step = p.angleincrement_rad;  % get data from struct
    
    th = ((0:1:size(RF,2)-1) * theta_step) + theta0;    % theta
    r = 1:1:size(RF,1);                                 % radius
    
    [im,~,~] = scanconvert(envelope,r,th);  % function given from professor
end


% Task 3b)
function showImage(im, figNumber, figTitle)
% Shows image
    dyn = 40;       % dynamic range (contrast)
    gain = -90;     % gain (brightness)
    
    figure(figNumber)
    imagesc(20*log10(abs(im)))
    clim([-dyn 0]-gain)
    colormap(gray)
    title(figTitle)
end


% Perform tasks
function IQ = performTask1(RF, f0, fs, figNumber, figTitle)
    Nfft = size(RF,1);
    freq = linspace(-0.5, 0.5, Nfft+1)*fs; freq(end) = [];  % create frequency array
    
    IQ = demodulation(RF, f0, fs);      % perform demodulation
    
    ps_RF = averagePower(RF);           % power spectra of RF
    ps_IQ = averagePower(IQ);           % power spectra of IQ
    
    
    figure(figNumber)
    plot(freq, ps_RF), hold on
    plot(freq, ps_IQ);
    xlabel('Frequency [Hz]'), ylabel('Amplitude [dB]')
    legend({'RF signal (before demodulation)'; 'IQ signal (after demodulation'})
    title(figTitle)
end

function envelope = performTask2(p, RF, fs, IQ, figNumber, figTitle)
    index = size(RF,2)/2;               % index of middle line
    [depth, line] = getLine(p, RF, fs, index);
    envelope = getEnvelope(IQ);
    
    
    figure(figNumber)
    plot(depth*1E3, line), hold on
    plot(depth*1E3, envelope(:,index), 'linewidth', 1.5)
    xlabel('Depth [mm]'), ylabel('Signal')
    legend({'Image line'; 'Envelope'})
    title(figTitle)
end

function performTask3(p, RF, envelope, figNumber, figTitle)
    im = beamToCartesian(p, RF, envelope);
    showImage(im, figNumber, figTitle);
end