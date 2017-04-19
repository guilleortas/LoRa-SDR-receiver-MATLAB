close all
clear

%% Created by Guillermo Ortas Delgado

%% Recorded signal import
[in, Fs] = audioread('lora-recorded-signal.wav');
% Allocate in-phase and quadrature components
x = complex(in(:,2), in(:,1)).';
clear in

%% LoRa parameters
BW = 125e3;
SF = 12;
Fc = 1.6385e6; % 1.64e6
symbol_time = 2^SF/BW; % 32.8e-3
symbols_per_frame = 57;

%% Signal channelization (DDC)
% Bring signal to baseband
t = 0:1/Fs:length(x)/Fs-1/Fs;
x = x.*cos(2*pi*Fc*t);
% Filter the signal
% freqs = [0 2*BW/Fs 2*BW/Fs*9/8 1];
% damps = [1 1 0 0];
% order = 50;
% b = firpm(order,freqs,damps);
% x = filter(b,1,x);
% Resampling
x = resample(x, BW, Fs); % Output sampling frequency is BW
Fs = BW;
clear t

%% Chirp generation
f0 = -BW/2;
f1 = BW/2;
t = 0:1/Fs:symbol_time - 1/Fs;

chirpI = chirp(t, f0, symbol_time, f1, 'linear', 90);
chirpQ = chirp(t, f0, symbol_time, f1, 'linear', 0);
upChirp = complex(chirpI, chirpQ);
clear chirpI chirpQ
% Normalization
% scale = max(max(abs(chirpI)), max(abs(chirpQ)));
% upChirp = upChirp / scale;

upChirp = repmat(upChirp,1,10);

%% Signal synchronization and cropping
% Find the start of the signal
[corr, lag] = xcorr(x, upChirp);
corrThresh = max(abs(corr))/4;
cLag = find(abs(corr) > corrThresh, 1);
signalStartIndex = abs(lag(cLag)) + 9*symbol_time*Fs;
signalEndIndex = round(signalStartIndex + symbols_per_frame*symbol_time*Fs);

% signalEndIndex = find(abs(corr) > corrThresh, length(corr));
% signalEndIndex = lag(signalEndIndex(end)) + 5*symbol_time*Fs;

% Synchronize SFD
symbol_offset = 0.25; % 12.25 to skip preamble and SFD
signalStartIndex = round(signalStartIndex + symbol_offset*symbol_time*Fs);
% Crop signal in time
% % signalStartIndex = 2.3055*Fs;
% % signalEndIndex = 4.11*Fs;
x = x(signalStartIndex:signalEndIndex);
clear lag corr

%% De-chirping
upChirp = repmat(upChirp,1,ceil(length(x)/length(upChirp)));
upChirp = upChirp(1:length(x));
de_chirped = x.*conj(upChirp);

%% Spectrogram computation
% To create a spectrogram 'grid' of symbols, these are the parameters.
signal = de_chirped;
Nfft = 2^SF; % 512
window_length = Nfft; % same as symbol_time*Fs;
[s, f, t] = spectrogram(signal, blackman(window_length), 0, Nfft, Fs);

%% Spectrogram plotting
surf(f,t,10*log10(abs(s.')),'EdgeColor','none')
axis xy; axis tight; colormap(jet); view(0,90);
ylabel('Time');
xlabel('Frequency (Hz)');
% xlim([0 BW])
ylim([0.0001 1])

%% Bit extraction
% s = s(:,1:symbols_per_frame-2);
[~, symbols] = max(abs(s));
symbols = mod(symbols - round(mean(symbols(1:8))), 2^SF);
% symbols = mod(symbols - 2^(SF-2) +1, 2^SF); % subject to that 0.25 symbol offset
bits =  dec2base(symbols, 2);
