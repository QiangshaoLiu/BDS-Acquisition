clc;
close all;

%data_in_11 = AD1_IN_IBUF90(1:end);   %Channel 1
%data_in_22 = AD2_IN_IBUF90(1:end);   %Channel 2

data_in_1 = csvread('L1_-40_12.csv',1,3,[1 3 131072 3]);%Channel 1
%data_in_2 = csvread('L1_-40_12.csv',1,4,[1 4 131072 4]);%Channel 2
data_in_1 = csvread('tianxian_14pm45_62.csv',1,3,[1 3 8192 3]);%Channel 2
[data1 data2 data3] = textread('IFtest.txt','%d%d%d');
data = [data1 data2 data3];

len = length(data3);
freq = (1:len)./len * 62;

figure(1)
plot(data3);

fft_data_squ = abs(fft(data3)).^2;
fft_data_log = 10*log10(fft_data_squ);

figure(2)
plot(freq,fft_data_log);   %带通采样，中频频点为9.22MHz和3.18MHz，带宽为4.092MHz

