%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Min Suk Lee
% Date: 8/16/2021
% Description: Load and format BioRadio csv file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
%% load in data from BR output file
%addpath('E:\Matlab_Git\BioRadio\HeadElectrode\NewNoise CSV Data')

% filenames %
filename31 = "Day19 Eyes Opened (Egg Electrodes).csv";
filename30 = "Day19 Eyes Opened.csv";
filename29 = "Day18 Eyes Opened (Egg Electrodes).csv";
filename28 = "Day18 Eyes Opened.csv";
filename27 = "Day17 Eyes Opened (Egg Electrodes).csv";
filename26 = "Day17 Eyes Opened.csv";
filename25 = "Day16 Eyes Opened (Egg Electrodes).csv";
filename24 = "Day16 Eyes Opened.csv";
filename23 = "Day15 Eyes Opened (Egg Electrodes).csv";
filename22 = "Day15 Eyes Opened.csv";
filename21 = "Day14 Eyes Opened v2 (Egg Electrodes).csv";
filename20 = "Day14 Eyes Opened.csv";
filename19 = "Day13 Eyes Opened (Egg Electrodes).csv";
filename18 = "Day13 Eyes Opened.csv";
filename17 = "Day12 Eyes Opened (Egg Electrodes).csv";
filename16 = "Day12 Eyes Opened.csv";
filename15 = "Day11 Eyes Opened (Egg Electrodes).csv";
filename14 = "Day11 Eyes Opened.csv";
filename13 = "Day10 Eyes Opened (Egg Electrodes).csv";
filename12 = "Day10 Eyes Opened.csv";
filename11 = "Day8 Eyes Opened (Egg Electrodes).csv";
filename10 = "Day8 Eyes Opened v2.csv";
filename9 = "Day7 Eyes Opened (Egg Electrodes).csv";
filename8 = "Day7 Eyes Opened.csv";
filename7 = "Day6 Eyes Opened v3 (Egg Electrodes).csv";
filename6 = "Day6 Eyes Opened v3.csv";
filename5 = "Day5 Eyes Opened.csv";
%filename4 = "Day4 Eyes Opened.csv";
filename3 = "Day3 Eyes Opened v2.csv";
filename2 = "Day2 Eyes Opened.csv";
filename1 = "Day1 Eyes Opened.csv";


%FILES = [filename1, filename2, filename3, filename4];
FILES = [filename1, filename2, filename3, filename5, filename6, filename7, ...
     filename8, filename9, filename10, filename11, filename12, filename13, ...
     filename14, filename15, filename16, filename17, filename18, filename19, ...
     filename20, filename21, filename22, filename23, filename24, filename25, ...
     filename26, filename27, filename28, filename29, filename30, filename31]; % Compile filenames %

LEGEND = FILES; % Create Filenames as Legends for graphing

% Extract data from file
for f = 1:length(FILES)
    file = readtable(FILES(f)); % get file as table
    if sum(sum(ismissing(file))) > 1 % check for any package drops
        disp("Package Drop: " + sum(sum(ismissing(file))));
    end
    NUMCHAN(f) = {width(file)-2}; % get number of channels
    EVENT_IDX(f) = {find(file.BioRadioEvent==1)}; % get event index
    DATA(f) = {file.(2)}; % get data
    TIMESTAMP(f) = {file.(1)}; % get time stamp
end
%% Initialize Variables %%
Fs = 1000; % Sampling Rate %

%% EEGlab Section
for f = 1:length(FILES)
    EEG(f) = pop_importdata('data', DATA{f}', 'nbchan', NUMCHAN{f}, 'srate', Fs);
end
%% Trim Data

% Marker Index for each file %
%str_mark = [1,1,1,1,1,1];
str_mark = ones(length(FILES),1);
%end_mark = [2,2];
str_offset = 5*Fs; % offset from the start index
duration = 30*Fs-1; % duration in seconds

for f = 1:length(FILES)
    str_idx = EVENT_IDX{f}(str_mark(f))+str_offset;
    sigtmp(f) = {DATA{f}(str_idx:(str_idx+duration))'};
    %time(f) = {TIMESTAMP{f}(str_idx:(str_idx+duration))'};
end
time = (1:(duration+1))/Fs;

%% Apply Highpass FIR Filter
hpFilt = designfilt('highpassfir','StopbandFrequency',0.25, ...
         'PassbandFrequency',0.35,'PassbandRipple',0.5, ...
         'StopbandAttenuation',65,'DesignMethod','kaiserwin');
%freqz(f1)
for f = 1:length(FILES)
    filterdata_1(f,:) = filtfilt(hpFilt,sigtmp{f});
end

%% Apply 60Hz Notch Filter
% design 60 Hz notch filter
f1 = designfilt('bandstopiir','FilterOrder',8, ...
               'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
               'DesignMethod','butter','SampleRate',Fs);
%freqz(f1)
for f = 1:length(FILES)
    filterdata(f,:) = filtfilt(f1,filterdata_1(f,:));
end 

%% Plot Time Series Data
figure(19);
hold on
for f = 1:length(FILES)
    scatter(time,filterdata(f,:).*10^6,'*','LineWidth',2);
end
box on
set(gca,'FontSize',14);
set(gca,'LineWidth',2);
xlabel('Time [sec]','FontWeight', 'bold','FontSize',15);
ylabel('Voltage [uV]','FontWeight', 'bold','FontSize',12);
lgd = legend(LEGEND,'Interpreter','none');
lgd.FontSize = 15;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
%ylim([-140 -60]);
%xlim([1 70]);

%% Plot Square Mean of Filtered Data
figure(18);
bar([1:length(FILES)],sqrt(mean((filterdata.*10^6).^2,2)));
%set(gca,'FontSize',14);
set(gca,'LineWidth',2);
set(gca,'XTick', 1:length(FILES));
set(gca,'xticklabel',FILES)
box on
xtickangle(45);
xlabel('File Type','FontWeight', 'bold','FontSize',15);
ylabel('Root Mean Square Voltage [uV]','FontWeight', 'bold','FontSize',12);


%% PSD
% Initialize PSD variables %
win_size = 2;
freqs_plt = [1 70];
toplot = 'OFF';

for f = 1:length(FILES)
    if length(str_mark)~=length(FILES)
        disp("str_mark length does not match the number of files")
        break;
    end
    %totsiz = EEG(f).pnts;
    [power, freqs] = spectopo(sigtmp{f}, 0, EEG(f).srate, ...
        'freqrange', freqs_plt, ...
        'winsize', Fs*win_size, ...
        'limits', [0 60 NaN NaN], ... %xmin xmax ymin ymax
        'plot', toplot);
    POWER(f) = {power};
    FREQ(f) = {freqs};
end

%% Plot PSD

figure(20);
hold on;
for f = 1:length(FILES)
    plot(FREQ{f},POWER{f},'LineWidth',2)
end
axis([0 70 -190 -130])
box on
set(gca,'FontSize',14);
set(gca,'LineWidth',2);
xlabel('Frequency [Hz]','FontWeight', 'bold','FontSize',15);
ylabel('PSD [10log_{10}(V^{2}/Hz)]','FontWeight', 'bold','FontSize',15);
lgd = legend(LEGEND,'Interpreter','none');
lgd.FontSize = 15;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
ylim([-140 -60]);
xlim([1 100]);