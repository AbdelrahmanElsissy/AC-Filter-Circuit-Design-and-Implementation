%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Exercise 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assignment 1
clc, close all; clear all;

% Initial Parameters-------------------------------------------------------
u1 = 1; %V
R = 2e3; %Ohms
C = 100e-9; %F
f = logspace(1,5,100);
f_c = 1 / (R*C*2*pi);
%--------------------------------------------------------------------------

% ---------------------Getting the Transfer Function-----------------------
% Magnitude----------------------------------------------------------------
u2 = (u1*i*2*pi*R*C*f) ./ (1 + i*2*pi*R*C*f);

figure;
subplot(2,1,1)
loglog(f/1000, abs(u2), 'LineWidth', 2);
xlabel("Frequency [kHz]"); ylabel("|H(f)|");
axis([20/1000 20 3e-2 2]);

H = line([f_c/1000 2.99],[1 1]);
set(H, 'LineStyle', '-.', 'Color', 'r');
H = line([f_c/1000 0.008],[1 0.01]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c/1000,1.3,'f_c', 'FontSize', 8, 'Color', 'r');
text(0.1,0.5,'-20 dB/Dek', 'FontSize', 8, 'Color', 'r');

% Phase--------------------------------------------------------------------
subplot(2,1,2)
semilogx(f/1000, angle(u2)*180/pi, 'LineWidth', 2);
xlabel("Frequency [kHz]"); ylabel("Phase (H(f)) [Deg]");
axis([20/1000 20 0 100]);
set(gca, 'YTick', [0 45 90])%, 'FontSize',14);
H= line([0.01 20], [45 45]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c/1000, 35, 'f_c', 'FontSize', 8, 'Color', 'r');
% -------------------------------------------------------------------------
%% Assignment 2
clc, clear all; close all;

% Initial Parameters-------------------------------------------------------
fs = 48000;
T = 4;
N = fs*T;
t = (0:N-1)'/fs;
sample_length = numel(t);
% -------------------------------------------------------------------------

% Signal Generation--------------------------------------------------------
signal = zeros(sample_length,1);
signal(sample_length/2) = 1; % Delta Peak of 1 a.u. at 2s.
%--------------------------------------------------------------------------

% Palying / Recording the Data---------------------------------------------
deviceReader = audioDeviceReader;
deviceReader.Driver = 'ASIO';
devices = getAudioDevices(deviceReader);

BufferSize = 1024;     % ASIO Buffer Size (Samples)
saving = false;

aPR = audioPlayerRecorder('Device', 'ASIO4ALL v2',...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);

[playData, recData, N_underrun,N_overrun] = play_rec(aPR, signal);
% -------------------------------------------------------------------------

t_recorded = (0:(numel(recData(:,1))-1))'/fs; % Time in (s) fo the Recording

% FFT of Channel 1---------------------------------------------------------
ch1_fft = fft(recData(:,1));
ch1_scale = abs(ch1_fft/numel(recData(:,1)));
ch1_positive = ch1_scale(1:numel(recData(:,1))/2+1);
ch1_positive(2:end-1) = 2*ch1_positive(2:end-1);

frequencies_ch1 = fs/numel(recData(:,1)) * (0:numel(recData(:,1))/2);
% -------------------------------------------------------------------------

% FFT of Channel 2---------------------------------------------------------
ch2_fft = fft(recData(:,2));
ch2_scale = abs(ch2_fft/numel(recData(:,2)));
ch2_positive = ch2_scale(1:numel(recData(:,2))/2+1);
ch2_positive(2:end-1) = 2*ch2_positive(2:end-1);

frequencies_ch2 = fs/numel(recData(:,2)) * (0:numel(recData(:,2))/2);
% -------------------------------------------------------------------------

% Plotting Signals (time-domain)-------------------------------------------
figure('Position', [100 100 400 400]);
plot(t_recorded, recData(:,1), '-b', 'Linewidth', 1); axis tight; grid on;
xlabel("Time (s)"); ylabel("Signal Amplitude (a.u.)");

figure('Position', [100 100 400 400]);
plot(t_recorded, recData(:,2), '-r', 'Linewidth', 1); axis tight; grid on;
xlabel("Time (s)"); ylabel("Signal Amplitude (a.u.)");
% -------------------------------------------------------------------------

% Ch2 of mic = red  
% Ch1 of mic = white

% ---------------------Getting the Transfer Function-----------------------
% Magnitude----------------------------------------------------------------
T_func = ch1_positive ./ ch2_positive ;
R = 2e3; %Ohms
C = 100e-9; %F
f_c = 1 / (R*C*2*pi);

figure ('Position', [100 100 800 400]);
subplot(2,1,1)
loglog(frequencies_ch1,T_func,"LineWidth", 1), axis tight; grid on;
xlabel("Frequency (Hz)"); ylabel("|H(f)|");
H = line([f_c f_c],[min(T_func) max(T_func)]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c+100 ,0.2,'f_c', 'FontSize', 8, 'Color', 'r');
% -------------------------------------------------------------------------

% Phase--------------------------------------------------------------------
subplot(2,1,2)

% Calulating the Raw FFT---------------------------------------------------
ch1_scale_raw = ch1_fft/numel(recData(:,1));
ch1_positive_raw = ch1_scale_raw(1:numel(recData(:,1))/2+1);
ch1_positive_raw(2:end-1) = 2*ch1_positive_raw(2:end-1);
ch2_scale_raw = ch2_fft/numel(recData(:,2));
ch2_positive_raw = ch2_scale_raw(1:numel(recData(:,2))/2+1);
ch2_positive_raw(2:end-1) = 2*ch2_positive_raw(2:end-1);
% -------------------------------------------------------------------------

angle = angle(ch1_positive_raw ./ ch2_positive_raw)*180/pi;
semilogx(frequencies_ch1, angle, 'LineWidth', 1); axis tight; grid on;
xlabel("Frequency [Hz]"); ylabel("Phase (H(f)) / Deg");
H = line([f_c f_c],[-200 200]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c+100 ,0.2,'f_c', 'FontSize', 8, 'Color', 'r');
% -------------------------------------------------------------------------
%% Assignment 3
clc, close all; clear all;

% Initial Parameters-------------------------------------------------------
fs = 48000;
T = 5;
N = fs*T;
t = (0:N-1)'/fs;
sample_length = numel(t);
hann_win = fs * 0.1; 
window = ramp_filter(hann_win, sample_length);
% -------------------------------------------------------------------------

% Generating the signal----------------------------------------------------
noise = wgn(sample_length, 1, 1,1,1,'linear') / 16;
signal= noise .* window * sample_length/sum(window); % White Noise
signal = signal(hann_win:sample_length-hann_win-1);
% -------------------------------------------------------------------------

% Palying / Recording the Data---------------------------------------------
deviceReader = audioDeviceReader;
deviceReader.Driver = 'ASIO';
devices = getAudioDevices(deviceReader);

BufferSize = 1024;     % ASIO Buffer Size (Samples)
saving = false;

aPR = audioPlayerRecorder('Device', 'ASIO4ALL v2',...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);

[playData, recData, N_underrun,N_overrun] = play_rec(aPR, signal);
% -------------------------------------------------------------------------

t_recorded = (0:(numel(recData(:,1))-1))'/fs; % Time in (s) fo the Recording

% FFT of Channel 1---------------------------------------------------------
ch1_fft = fft(recData(:,1));
ch1_scale = abs(ch1_fft/numel(recData(:,1)));
ch1_positive = ch1_scale(1:numel(recData(:,1))/2+1);
ch1_positive(2:end-1) = 2*ch1_positive(2:end-1);

frequencies_ch1 = fs/numel(recData(:,1)) * (0:numel(recData(:,1))/2);
% -------------------------------------------------------------------------

% FFT of Channel 2---------------------------------------------------------
ch2_fft = fft(recData(:,2));
ch2_scale = abs(ch2_fft/numel(recData(:,2)));
ch2_positive = ch2_scale(1:numel(recData(:,2))/2+1);
ch2_positive(2:end-1) = 2*ch2_positive(2:end-1);

frequencies_ch2 = fs/numel(recData(:,2)) * (0:numel(recData(:,2))/2);
% -------------------------------------------------------------------------

% Plotting Signals (time-domain)-------------------------------------------
figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,1), '-b', 'Linewidth', 1); axis tight, grid on;
xlabel("Time (s)"); ylabel("Signal Amplitude (a.u.)");

figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,2), '-r', 'Linewidth', 1); axis tight, grid on;
xlabel("Time (s)"); ylabel("Signal Amplitude (a.u.)");
% -------------------------------------------------------------------------

% Ch2 of mic = red  
% Ch1 of mic = white

% ------------------Getting the Trasfer Function---------------------------
% Magnitude----------------------------------------------------------------
T_func = ch1_positive ./ ch2_positive ;
R = 2e3; %Ohms
C = 100e-9; %F
f_c = 1 / (R*C*2*pi);

figure('Position', [100 100 800 400]);
subplot(2,1,1)
loglog(frequencies_ch1,T_func,"LineWidth", 1), axis tight; grid on;
xlabel("Frequency [Hz]"); ylabel("|H(f)|");
H = line([f_c f_c],[min(T_func) max(T_func)]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c+100 ,0.2,'f_c', 'FontSize', 8, 'Color', 'r');
% -------------------------------------------------------------------------

% Phase--------------------------------------------------------------------
subplot(2,1,2)

% Calulating the Raw FFT---------------------------------------------------
ch1_scale_raw = ch1_fft/numel(recData(:,1));
ch1_positive_raw = ch1_scale_raw(1:numel(recData(:,1))/2+1);
ch1_positive_raw(2:end-1) = 2*ch1_positive_raw(2:end-1);
ch2_scale_raw = ch2_fft/numel(recData(:,2));
ch2_positive_raw = ch2_scale_raw(1:numel(recData(:,2))/2+1);
ch2_positive_raw(2:end-1) = 2*ch2_positive_raw(2:end-1);
% -------------------------------------------------------------------------

angle = angle(ch1_positive_raw ./ ch2_positive_raw)*180/pi;
semilogx(frequencies_ch1, angle, 'LineWidth', 1); axis tight; grid on;
xlabel("Frequency [Hz]"); ylabel("Phase (H(f)) / Deg");
H = line([f_c f_c],[-200 200]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c+100 ,0.2,'f_c', 'FontSize', 8, 'Color', 'r');
% -------------------------------------------------------------------------
%% Assignment 4
clc, clear all; close all;

% Initial Parameters-------------------------------------------------------
startFreq = 100;
endFreq = 24000;
T = 10;
fs = 48000;
BufferSize = 1024;
t = 0:1/fs:T;
% -------------------------------------------------------------------------

% Generating the signal----------------------------------------------------
signal = chirp(t, startFreq, T, endFreq, 'linear');
signal = signal * 0.2; % Linear-sweep Signal
% -------------------------------------------------------------------------

% Palying / Recording the Data---------------------------------------------
deviceReader = audioDeviceReader;
deviceReader.Driver = 'ASIO';
devices = getAudioDevices(deviceReader);


saving = false;

aPR = audioPlayerRecorder('Device', 'ASIO4ALL v2',...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);

[playData, recData, N_underrun,N_overrun] = play_rec(aPR, signal');
% -------------------------------------------------------------------------

t_recorded = (0:(numel(recData(:,1))-1))'/fs; % Time in (s) fo the Recording

% FFT of Channel 1---------------------------------------------------------
ch1_fft = fft(recData(:,1));
ch1_scale = abs(ch1_fft/numel(recData(:,1)));
ch1_positive = ch1_scale(1:numel(recData(:,1))/2+1);
ch1_positive(2:end-1) = 2*ch1_positive(2:end-1);

frequencies_ch1 = fs/numel(recData(:,1)) * (0:numel(recData(:,1))/2);
% -------------------------------------------------------------------------

% FFT of Channel 2---------------------------------------------------------
ch2_fft = fft(recData(:,2));
ch2_scale = abs(ch2_fft/numel(recData(:,2)));
ch2_positive = ch2_scale(1:numel(recData(:,2))/2+1);
ch2_positive(2:end-1) = 2*ch2_positive(2:end-1);

frequencies_ch2 = fs/numel(recData(:,2)) * (0:numel(recData(:,2))/2);
% -------------------------------------------------------------------------

% Plotting Signals (time-domain)-------------------------------------------
figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,1), '-b', 'Linewidth', 1); grid on, axis tight;
xlabel("Time (s)"); ylabel("Signal Amplitude (a.u.)");

figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,2), '-r', 'Linewidth', 1); grid on, axis tight;
xlabel("Time (s)"); ylabel("Signal Amplitude (a.u.)");
% -------------------------------------------------------------------------

% Plotting FFT's-----------------------------------------------------------
figure('Position', [100 100 800 400]);
loglog(frequencies_ch1,ch1_positive, "LineWidth", 1), axis tight; grid on;
xlabel("Frequency (Hz)");
ylabel("Signal Amplitude (a.u.)");

figure('Position',[100 100 800 400]);
loglog(frequencies_ch2,ch2_positive, "LineWidth", 1), axis tight; grid on;
xlabel("Frequency (Hz)");
ylabel("Signal Amplitude (a.u.)");
% -------------------------------------------------------------------------

% Ch2 of mic = red  
% Ch1 of mic = white

% ------------------Getting the Trasfer Function---------------------------
% Magnitude----------------------------------------------------------------
T_func = ch1_positive ./ ch2_positive ;
R = 2e3; %Ohms
C = 100e-9; %F
f_c = 1 / (R*C*2*pi);

figure('Position',[100 100 800 400]);
subplot(2,1,1)
loglog(frequencies_ch2,T_func,"LineWidth", 1), axis tight; grid on;
xlabel("Frequency [Hz]"); ylabel("|H(f)|");
H = line([f_c f_c],[min(T_func) max(T_func)]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c+100 ,0.2,'f_c', 'FontSize', 8, 'Color', 'r');
% -------------------------------------------------------------------------

% Phase--------------------------------------------------------------------
subplot(2,1,2)

% Calulating the Raw FFT---------------------------------------------------
ch1_scale_raw = ch1_fft/numel(recData(:,1));
ch1_positive_raw = ch1_scale_raw(1:numel(recData(:,1))/2+1);
ch1_positive_raw(2:end-1) = 2*ch1_positive_raw(2:end-1);
ch2_scale_raw = ch2_fft/numel(recData(:,2));
ch2_positive_raw = ch2_scale_raw(1:numel(recData(:,2))/2+1);
ch2_positive_raw(2:end-1) = 2*ch2_positive_raw(2:end-1);
% -------------------------------------------------------------------------

angle = angle(ch1_positive_raw ./ ch2_positive_raw)*180/pi;
semilogx(frequencies_ch2, angle, 'LineWidth', 1); axis tight; grid on;
xlabel("Frequency [Hz]"); ylabel("Phase (H(f)) / Deg");
H = line([f_c f_c],[-200 200]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c+100 ,0.2,'f_c', 'FontSize', 8, 'Color', 'r');
% -------------------------------------------------------------------------
%% Assignment 5
clc, clear all, close all;

% Initial Parameters-------------------------------------------------------
startFreq = 100;
endFreq = 24000;
T = 10;
fs = 48000;
BufferSize = 1024;
t = 0:1/fs:T;
% -------------------------------------------------------------------------

% Generating the signal----------------------------------------------------
signal = chirp(t, startFreq, T, endFreq, 'logarithmic');
signal = signal * 0.2; % Log-sweep Signal
% -------------------------------------------------------------------------

% Palying / Recording the Data---------------------------------------------
deviceReader = audioDeviceReader;
deviceReader.Driver = 'ASIO';
devices = getAudioDevices(deviceReader);


saving = false;

aPR = audioPlayerRecorder('Device', 'ASIO4ALL v2',...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);

[playData, recData, N_underrun,N_overrun] = play_rec(aPR, signal');
% -------------------------------------------------------------------------

t_recorded = (0:(numel(recData(:,1))-1))'/fs; % Time in (s) fo the Recording

% FFT of Channel 1---------------------------------------------------------
ch1_fft = fft(recData(:,1));
ch1_scale = abs(ch1_fft/numel(recData(:,1)));
ch1_positive = ch1_scale(1:numel(recData(:,1))/2+1);
ch1_positive(2:end-1) = 2*ch1_positive(2:end-1);

frequencies_ch1 = fs/numel(recData(:,1)) * (0:numel(recData(:,1))/2);
% -------------------------------------------------------------------------

% FFT of Channel 2---------------------------------------------------------
ch2_fft = fft(recData(:,2));
ch2_scale = abs(ch2_fft/numel(recData(:,2)));
ch2_positive = ch2_scale(1:numel(recData(:,2))/2+1);
ch2_positive(2:end-1) = 2*ch2_positive(2:end-1);

frequencies_ch2 = fs/numel(recData(:,2)) * (0:numel(recData(:,2))/2);
% -------------------------------------------------------------------------

% Plotting Signals (time-domain)-------------------------------------------
figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,1), '-b', 'Linewidth', 1); grid on, axis tight;
xlabel("Time (s)"); ylabel("Signal Amplitude (a.u.)");

figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,2), '-r', 'Linewidth', 1); grid on, axis tight;
xlabel("Time (s)"); ylabel("Signal Amplitude (a.u.)");
% -------------------------------------------------------------------------

% Plotting FFT's-----------------------------------------------------------
figure('Position', [100 100 800 400]);
loglog(frequencies_ch1,ch1_positive, "LineWidth", 1), axis tight; grid on;
xlabel("Frequency (Hz)");
ylabel("Signal Amplitude (a.u.)");

figure('Position',[100 100 800 400]);
loglog(frequencies_ch2,ch2_positive, "LineWidth", 1), axis tight; grid on;
xlabel("Frequency (Hz)");
ylabel("Signal Amplitude (a.u.)");
% -------------------------------------------------------------------------

% Ch2 of mic = red  
% Ch1 of mic = white

% ------------------Getting the Trasfer Function---------------------------
% Magnitude----------------------------------------------------------------
T_func = ch1_positive ./ ch2_positive ;
R = 2e3; %Ohms
C = 100e-9; %F
f_c = 1 / (R*C*2*pi);

figure('Position', [100 100 800 400]);
subplot(2,1,1)
loglog(frequencies_ch2,T_func,"LineWidth", 1), axis tight; grid on;
title("Signal in Frequency-Domain (Logaritmic Scale)");
xlabel("Frequency [Hz]"); ylabel("|H(f)|");
H = line([f_c f_c],[min(T_func) max(T_func)]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c+100 ,0.2,'f_c', 'FontSize', 8, 'Color', 'r');

% Phase--------------------------------------------------------------------
subplot(2,1,2)

% Calulating the Raw FFT---------------------------------------------------
ch1_scale_raw = ch1_fft/numel(recData(:,1));
ch1_positive_raw = ch1_scale_raw(1:numel(recData(:,1))/2+1);
ch1_positive_raw(2:end-1) = 2*ch1_positive_raw(2:end-1);
ch2_scale_raw = ch2_fft/numel(recData(:,2));
ch2_positive_raw = ch2_scale_raw(1:numel(recData(:,2))/2+1);
ch2_positive_raw(2:end-1) = 2*ch2_positive_raw(2:end-1);
% -------------------------------------------------------------------------

angle = angle(ch1_positive_raw ./ ch2_positive_raw)*180/pi;
semilogx(frequencies_ch2, angle, 'LineWidth', 1); axis tight; grid on;
xlabel("Frequency [Hz]"); ylabel("Phase (H(f)) / Deg");
H = line([f_c f_c],[-200 200]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c+100 ,0.2,'f_c', 'FontSize', 8, 'Color', 'r');
% -------------------------------------------------------------------------
%% Assingment 6
clc, close all, clear all;

% Initial Parameters-------------------------------------------------------
fs = 48000;
T = 1;
N = fs*T;
t = (0:N-1)'/fs;
BufferSize = 1024;
sample_length = numel(t);
hann_win = fs * 0.1; 
window = ramp_filter(hann_win, sample_length);
frequencies = [10 30 100 300 1e3 3e3 10e3 20e3];
% -------------------------------------------------------------------------

% Generating the signal----------------------------------------------------
for j = 1:length(frequencies)
    signals(:, j) = 0.1*cos((2*pi)*t*frequencies(j));
    ramp_sig(:, j) = signals(:, j) .* window * sample_length/sum(window);
end
% -------------------------------------------------------------------------

% Selecting a Portion of the Signal----------------------------------------
signal1 = signals(:,1);
signal1 = signal1(hann_win:sample_length-hann_win-1);
signal2 = signals(:,2);
signal2 = signal2(hann_win:sample_length-hann_win-1);
signal3 = signals(:,3);
signal3 = signal3(hann_win:sample_length-hann_win-1);
signal4 = signals(:,4);
signal4 = signal4(hann_win:sample_length-hann_win-1);
signal5 = signals(:,5);
signal5 = signal5(hann_win:sample_length-hann_win-1);
signal6 = signals(:,6);
signal6 = signal6(hann_win:sample_length-hann_win-1);
signal7 = signals(:,7);
signal7 = signal7(hann_win:sample_length-hann_win-1);
signal8 = signals(:,8);
signal8 = signal8(hann_win:sample_length-hann_win-1);
% -------------------------------------------------------------------------

% Palying / Recording the Data---------------------------------------------
deviceReader = audioDeviceReader;
deviceReader.Driver = 'ASIO';
devices = getAudioDevices(deviceReader);


saving = false;

aPR = audioPlayerRecorder('Device', 'ASIO4ALL v2',...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);


[playData1, recData1, N_underrun1, N_overrun1] = play_rec(aPR, signal1);
t1_recorded = (0:(numel(recData1(:,1))-1))'/fs;

[playData2, recData2, N_underrun2, N_overrun2] = play_rec(aPR, signal2);
t2_recorded = (0:(numel(recData2(:,1))-1))'/fs;

[playData3, recData3, N_underrun3, N_overrun3] = play_rec(aPR, signal3);
t3_recorded = (0:(numel(recData3(:,1))-1))'/fs;

[playData4, recData4, N_underrun4, N_overrun4] = play_rec(aPR, signal4);
t4_recorded = (0:(numel(recData4(:,1))-1))'/fs;

[playData5, recData5, N_underrun5, N_overrun5] = play_rec(aPR, signal5);
t5_recorded = (0:(numel(recData1(:,1))-1))'/fs;

[playData6, recData6, N_underrun6, N_overrun6] = play_rec(aPR, signal6);
t6_recorded = (0:(numel(recData1(:,1))-1))'/fs;

[playData7, recData7, N_underrun7, N_overrun7] = play_rec(aPR, signal7);
t7_recorded = (0:(numel(recData7(:,1))-1))'/fs;

[playData8, recData8, N_underrun8, N_overrun8] = play_rec(aPR, signal8);
t8_recorded = (0:(numel(recData8(:,1))-1))'/fs;
% -------------------------------------------------------------------------

% Grouping Ch1-------------------------------------------------------------
rec_data_ch1 = [recData1(:,1),recData2(:,1),recData3(:,1),...
                recData4(:,1),recData5(:,1),recData6(:,1),...
                recData7(:,1),recData8(:,1)];
% -------------------------------------------------------------------------

% Grouping Ch2-------------------------------------------------------------
rec_data_ch2 = [recData1(:,2),recData2(:,2),recData3(:,2),...
                recData4(:,2),recData5(:,2),recData6(:,2),...
                recData7(:,2),recData8(:,2)];
% -------------------------------------------------------------------------

% Grouping T_recorded------------------------------------------------------
time_rec = [t1_recorded,t2_recorded,t3_recorded,t4_recorded,t5_recorded,...
            t6_recorded,t7_recorded,t8_recorded];
% -------------------------------------------------------------------------

% FFT of Channel 1---------------------------------------------------------
for m = 1:8
        signal = rec_data_ch1(:,m);

        signal_length = numel(signal);
      
        fft_signal = fft(signal);
        fft_scaled = fft_signal/signal_length;

        fft_magnitude = abs(fft_scaled);
        fft_positive = fft_magnitude(1:signal_length/2+1);
        fft_positive(2:end-1) = 2*fft_positive(2:end-1);
        
        fft_positive_raw = fft_scaled(1:signal_length/2+1);
        fft_positive_raw(2:end-1) = 2*fft_positive_raw(2:end-1);

        frequencies = fs/signal_length * (0:(signal_length/2));

        multiple_fft_ch1(:, m) = fft_positive;
        multiple_fft_ch1_raw(:, m) = fft_positive_raw;
        multiple_frequencies_ch1(:, m) = frequencies;
end
% -------------------------------------------------------------------------

% FFT of Channel 2---------------------------------------------------------
for m = 1:8
        signal = rec_data_ch2(:,m);

        signal_length = numel(signal);
      
        fft_signal = fft(signal);
        fft_scaled = fft_signal/signal_length;

        fft_magnitude = abs(fft_scaled);
        fft_positive = fft_magnitude(1:signal_length/2+1);
        fft_positive(2:end-1) = 2*fft_positive(2:end-1);

        fft_positive_raw = fft_scaled(1:signal_length/2+1);
        fft_positive_raw(2:end-1) = 2*fft_positive_raw(2:end-1);
        
        frequencies = fs/signal_length * (0:(signal_length/2));

        multiple_fft_ch2(:, m) = fft_positive;
        multiple_frequencies_ch2(:, m) = frequencies;
        multiple_fft_ch2_raw(:, m) = fft_positive_raw;

        angles(:,m) = angle(multiple_fft_ch1_raw(:, m) ./ multiple_fft_ch2_raw(:, m))*180/pi;
end
% -------------------------------------------------------------------------

% Plotting Signals (time-domain)-------------------------------------------

fig1 = figure('Position', [20 20 800 800]);
% xlabel("Time (s)"); ylabel("Signal Amplitude (a.u.)");
for n = 1:8
    subplot(8,1,n)
    plot(time_rec(:,n), rec_data_ch1(:,n), '-b', 'Linewidth', 1); grid on, axis tight;   
end
han=axes(fig1,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel("Time (s)", 'FontSize',14); ylabel("Signal Amplitude (a.u.)", 'FontSize',14);

fig2 = figure('Position', [20 20 800 800]);
for n = 1:8
    subplot(8,1,n)
    plot(time_rec(:,n), rec_data_ch2(:,n), '-r', 'Linewidth', 1); grid on, axis tight;
end
han=axes(fig2,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel("Time (s)", 'FontSize',14); ylabel("Signal Amplitude (a.u.)", 'FontSize',14);
% -------------------------------------------------------------------------

% ------------------Getting the Trasfer Function---------------------------
R = 2e3; %Ohms
C = 100e-9; %F
f_c = 1 / (R*C*2*pi); % Cut-off Frequency
% Magnitude and Phase------------------------------------------------------
for m = 1:8
    T_func(:,m) = multiple_fft_ch1(:,m) ./ multiple_frequencies_ch2(:,m);

    figure;
    subplot(2,1,1)
    loglog(multiple_frequencies_ch2(:,m),T_func(:,m),"LineWidth", 1), axis tight; grid on;
    xlabel("Frequency [Hz]"); ylabel("|H(f)|");
    H = line([f_c f_c],[min(T_func(:,m)) 8e-5]);
    set(H, 'LineStyle', '-.', 'Color', 'r');

    subplot(2,1,2)
    semilogx(multiple_frequencies_ch2(:,m), angles(:,m), 'LineWidth', 1); axis tight; grid on;
    xlabel("Frequency [Hz]"); ylabel("Phase (H(f)) / Deg");
    H = line([f_c f_c],[-200 200]);
    set(H, 'LineStyle', '-.', 'Color', 'r');
end
% -------------------------------------------------------------------------
%% Assignment 7
clc, close all, clear all;

% Initial Parameters-------------------------------------------------------
T = 3.5;  % Total signal duration in (s), more averages with longer signal
t_start = 0.4;                   % signal start: remove also transient part
t_stop=3.4;                      % signal end: remove the fading part
% -------------------------------------------------------------------------

% Define parameters of your measurement------------------------------------
name='Test Transfer Function Measurement V0';  % name for output file
fs = 48000;             % sampling frequency, make sure value is correct!
N_freqs = 179;          % number of frequencies between f_start and f_stop 
f_start = 55;           % start frequency (Hz)
f_stop = 22000;         % stop frequency (Hz)
block_size = 2^12;      % block size (FFT window length)
ampl=0.01;              % select peak amplitude of output re full scale
% -------------------------------------------------------------------------

% Initialize audio interface and device reader-----------------------------
deviceReader = audioDeviceReader;
deviceReader.Driver = 'ASIO';

BufferSize = 1024;
saving = false;
aPR = audioPlayerRecorder('Device', 'ASIO4ALL v2', ...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);
%--------------------------------------------------------------------------

% The number of output frequencies might be less than N_freqs, as the function removes redundant frequencies  
[t, sig, ms_indices] = generate_multisine(N_freqs, f_start, f_stop, fs, block_size, T);
sig=sig*ampl;  % scale signal

% Ramp-up, ramp-down with hanning window-----------------------------------
t_ramp = 200e-3;
n_ramp = floor(t_ramp*fs);
hann = hanning(2*n_ramp, 'periodic');
sig(1:n_ramp) = sig(1:n_ramp) .* hann(1:n_ramp);
sig(end-n_ramp+1:end) = sig(end-n_ramp+1:end) .* hann(n_ramp+1:end);
% -------------------------------------------------------------------------

% Palying / Recording the Data---------------------------------------------
[playData, recData, N_underrun, N_overrun] = play_rec(aPR, sig);
release(aPR);

t_recData = (0:size(recData, 1)-1)'/fs;
% -------------------------------------------------------------------------

% Process measurement data-------------------------------------------------
n_start = floor(t_start*fs);
n_stop = floor(t_stop*fs);
num_avg = floor((n_stop-n_start+1) /block_size); % number of blocks in data
n_stop = n_start + num_avg*block_size;          % align end to block border

rec = recData(n_start:n_stop-1,:);  % cut out data for analysis
t_rec = (0:length(rec)-1)/fs;

fprintf('Analyse data from T=%.2f to %.2f (averaging over %i blocks)\n',t_start,t_stop,num_avg);
% -------------------------------------------------------------------------

% Average in time domain---------------------------------------------------
rec_avg=mean(reshape(rec,block_size,num_avg,2),2);
% -------------------------------------------------------------------------

% FFT of input and output--------------------------------------------------
fft_rec = fft(rec_avg)/block_size;
fft_rec = 2*fft_rec(2:floor(block_size/2)+1,:); % single sided spectrum without DC
fft_freqs = fs*(1:block_size/2)'/block_size; % frequencies of spectrum
% -------------------------------------------------------------------------
    
% Spectrum for frequencies, where energy was provided----------------------
ms_fft_freqs = fft_freqs(ms_indices,:); % select frequency vector
ms_fft_rec = fft_rec(ms_indices,:);     % select frequencies in spectrum
% -------------------------------------------------------------------------

% Calculate transfer function----------------------------------------------
H=ms_fft_rec(:, 1) ./ ms_fft_rec(:, 2);  % Modified due to our connections

R = 2e3; %Ohms
C = 100e-9; %F
f_c = 1 / (R*C*2*pi); % Cut-off Frequency

% Plotting the Transfer Fucntion-------------------------------------------
figure('Position',[100 100 800 400]);
subplot(2,1,1)
title('Transfer Function');
loglog(ms_fft_freqs, abs(H), 'b.');
L_fc = line([f_c f_c],[0.04 1]);
set(L_fc, 'LineStyle', '-.', 'Color', 'r');
text(f_c+50 ,0.3,'f_c', 'FontSize', 8, 'Color', 'r');
ylabel('|H|');

subplot(2,1,2)
semilogx(ms_fft_freqs, unwrap(angle(H))*180/pi, 'b.');
xlabel('f / Hz');
ylabel('Phase(H) / Deg');
L_fc_an= line([57 22000], [45 45]);
set(L_fc_an, 'LineStyle', '-.', 'Color', 'r');
text(f_c+800 , 32, 'f_c', 'FontSize', 8, 'Color', 'r');
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%