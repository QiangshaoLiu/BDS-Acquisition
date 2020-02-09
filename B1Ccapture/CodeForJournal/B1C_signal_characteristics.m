%本程序给出了北斗B1C信号的功率谱密度、BOC以及QMBOC的自相关函数仿真图
%作者：LSQ
%日期：2019年3月20日

clc
clear all
close all

load('C.mat');
load('W.mat');

f_medium = 15.58e6;      %中频频率
f_sample = 30*1.023e6;      %采样频率
f_code = 1.023e6;           %码速率
f0 = 1.023e6;               %基准频率
f1 = f0;                    %BOC(1,1)
f2 = 6*f0;                  %BOC(6,1)
sim_t = 1e-3;               %仿真时间
t = 0 : 1/f_sample : sim_t - 1/f_sample;
t_n = f_code * sim_t;       %仿真时间内的伪码数量
c_fs = f_sample / f_code;   %每个码元所被采样次数

for k = 1 : t_n
    i = mod(k, 18414000);
    if i == 0
        PseduCode(1+(k-1)*f_sample/f_code : k*f_sample/f_code) = C(18414000);
    else
        PseduCode(1+(k-1)*f_sample/f_code : k*f_sample/f_code) = C(i);
    end
end

SC1 = sign(sin(2*pi*f1*t));
SC2 = sign(sin(2*pi*f2*t));

S1 = PseduCode.*SC1.*cos(2*pi*f_medium*t);
S2 = PseduCode.*SC2.*cos(2*pi*f_medium*t);
S3 = PseduCode.*cos(2*pi*f_medium*t);

S1_fft = fft(S1);
S1_pow = abs(S1_fft).^2;
S2_fft = fft(S2);
S2_pow = abs(S2_fft).^2;
S3_fft = fft(S3);
S3_pow = abs(S3_fft).^2;

%%BOC信号的功率谱密度分析
figure(1)
plot(10*log10(S1_pow),'b');
hold on;
plot(10*log10(S2_pow),'g');
hold on;
plot(10*log10(S3_pow),'r');

xlabel('Frequency(kHz)');
ylabel('PSD(dBw)');
legend('BOC(1,1)','BOC(6,1)','BPSK');
%title('BOC信号特性分析');
grid on;
axis([0 3e4 -20 70]);

t_corr = 0 : 1/f_sample : 1/f_code;   %自相关函数取值范围为-1至1码片，需要一个码片的数据做自相关运算
SC1_prn = sign(sin(2*pi*f1*t_corr));  %BOC调制信号的自相关函数与扩频符号的自相关函数相等
SC2_prn = sign(sin(2*pi*f2*t_corr));
t_corr_length = length(t_corr);
Bpsk_one = ones(1,floor(t_corr_length));
%Bpsk_zero = zeros(1,floor(t_corr_length/2));
SC3_prn = [0 Bpsk_one];

[SC1_corr lag_1] = xcorr(SC1_prn, 'coeff');
[SC2_corr lag_2] = xcorr(SC2_prn, 'coeff');
[SC3_corr lag_3] = xcorr(SC3_prn, 'coeff');

index_1 = 2*lag_1/length(lag_1);
index_2 = 2*lag_2/length(lag_2);
index_3 = 2*lag_3/length(lag_2);

%%MBOC信号自相关性能分析
SC_mboc_corr = 10/11*SC1_corr + 1/11*SC2_corr;

%%QMBOC信号
SC_Qmboc_corr = 29/33*SC1_corr + 4/33*SC2_corr;  %BDS B1C信号

SB1C = 1/2*PseduCode.*SC1.*cos(2*pi*f_medium*t)+...
    sqrt(1/11)*PseduCode.*SC2.*cos(2*pi*f_medium*t)+...
    sqrt(-1)*sqrt(29/44)*PseduCode.*SC1.*cos(2*pi*f_medium*t);

SB1C_fft = fft(SB1C);
SB1C_pow = abs(SB1C_fft).^2;

alfa = 1;
alfa2 = 6;
beta = 1;
beta2 = 1;
beta_bpsk = 2;
pi = 3.1415926535898;
j = sqrt(-1);
fs = alfa*f0;
fs2 = alfa2*f0;
n = 2*alfa/beta;         %BOC(1,1)调制阶数
n2 = 2*alfa2/beta2;      %BOC(6,1)调制阶数
ts = 1/(2*fs);
ts2 = 1/(2*fs2);
tc = 1./(beta.*f0);
tc2 = 1./(beta2.*f0);
tc_bpsk = 1./(beta_bpsk.*f0);
f = -20./tc:0.001./tc:20./tc;
f2 = -20./tc2:0.001./tc2:20./tc2;
f_bpsk = -20./tc_bpsk:0.001./tc_bpsk:20./tc_bpsk;

f_x_boc = n.*ts.*sinc(f.*ts).^2.*sin(pi.*f.*ts).^2;
f_x_bpsk = tc_bpsk.*sinc(f_bpsk.*tc_bpsk).^2;
f_x_boc2 = 16.*ts2./n2.*sinc(f2.*ts2).^2.*sin(3.*pi.*f2.*ts2).^2.*(4.*cos(pi.*f2.*ts2).^2-3).^2.*cos(6.*pi.*f2.*ts2).^2;
f_mboc = 29./33.*f_x_boc + 4./33.*f_x_boc2;
f_qmboc = 29./33.*f_x_boc + j.*4./33.*f_x_boc2;

figure(2)
subplot(2,2,1);
plot(index_3, SC3_corr, 'g',index_1, SC1_corr,'b', index_2, SC2_corr,'r');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1 1 -0.6 1]);
legend('BPSK-R(2)','BOC(1,1)','BOC(6,1)');
subplot(2,2,2);
plot(index_1, SC_Qmboc_corr, index_2, SC1_corr,'r');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1 1 -0.6 1]);
legend('QMBOC(6,1,4/33)','BOC(1,1)');
subplot(2,2,3);
plot(f_bpsk,10.*log10(f_x_bpsk),'g-');axis([-8./tc_bpsk,8./tc_bpsk,-110,-60]);
hold on
plot(f,10.*log10(f_x_boc),'b-');axis([-8./tc,8./tc,-110,-60]);
hold on
plot(f2,10.*log10(f_x_boc2),'r-');axis([-8./tc2,8./tc2,-110,-60]);
xlabel('频率  Hz');
ylabel('功率谱密度 dBW/Hz');
legend('BPSK-R(2)','BOC(1,1)','BOC(6,1)');
subplot(2,2,4);
plot(f,10.*log10(f_x_boc),'r-');axis([-20./tc,20./tc,-110,-60]);
hold on
plot(f2,10.*log10(f_qmboc),'b-');axis([-20./tc2,20./tc2,-110,-60]);
xlabel('频率  Hz');
ylabel('功率谱密度 dBW/Hz');
legend('BOC(1,1)','QMBOC(6,1,4/33)');

figure(3)
plot(f,10.*log10(f_x_boc),'r-');axis([-20./tc,20./tc,-110,-60]);
hold on
plot(f2,10.*log10(f_qmboc),'b-');axis([-20./tc2,20./tc2,-110,-60]);
xlabel('Frequency(Hz)');
ylabel('PSD(dBW/Hz)');
legend('BOC(1,1)','QMBOC(6,1,4/33)');
grid on;

figure(4)
plot(index_1, SC_Qmboc_corr,'r', index_2, SC1_corr,index_2, SC2_corr, 'g');
xlabel('Code Delay(Chips)');
ylabel('Correlation Function');
axis([-1 1 -0.6 1]);
legend('QMBOC(6,1,4/33)','BOC(1,1)','BOC(6,1)');
grid on;








