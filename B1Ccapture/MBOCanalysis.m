%�ó�����Ҫ������MBOC�źŵ���������Ժ͹������ܶ�
clc
clear all
close all

load('C.mat');
load('W.mat');

%f_medium = 20*1.023e6;      %��ƵƵ��
%f_sample = 70*1.023e6;      %����Ƶ��
f_medium = 15.58e6;      %��ƵƵ��
f_sample = 70*1.023e6;      %����Ƶ��
f_code = 1.023e6;           %������
f0 = 1.023e6;               %��׼Ƶ��
f1 = f0;                    %BOC(1,1)
f2 = 6*f0;                  %BOC(6,1)
SNR = 20;                   %�����20db
sim_t = 1e-3;               %����ʱ��
t = 0 : 1/f_sample : sim_t - 1/f_sample;
t_n = f_code * sim_t;       %����ʱ���ڵ�α������
c_fs = f_sample / f_code;   %ÿ����Ԫ������������

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

%S1_add = awgn(S1, SNR);
%S2_add = awgn(S2, SNR);
%S3_add = awgn(S3, SNR);

S1_fft = fft(S1);
S1_pow = abs(S1_fft).^2;
S2_fft = fft(S2);
S2_pow = abs(S2_fft).^2;
S3_fft = fft(S3);
S3_pow = abs(S3_fft).^2;

%%BOC�źŵĹ������ܶȷ���
figure(1)
plot(10*log10(S1_pow),'b');
hold on;
plot(10*log10(S2_pow),'g');
hold on;
plot(10*log10(S3_pow),'r');

xlabel('Frequency(kHz)');
ylabel('PSD(dBw)');
legend('BOC(1,1)','BOC(6,1)','BPSK');
title('BOC�ź����Է���');

%%BOC�ź���������ܷ���
t_corr = 0 : 1/f_sample : 1/f_code;   %����غ���ȡֵ��ΧΪ-1��1��Ƭ����Ҫһ����Ƭ�����������������
SC1_prn = sign(sin(2*pi*f1*t_corr));  %BOC�����źŵ�����غ�������Ƶ���ŵ�����غ������
SC2_prn = sign(sin(2*pi*f2*t_corr));
t_corr_length = length(t_corr);
Bpsk_one = ones(1,floor(t_corr_length/2));
Bpsk_zero = zeros(1,floor(t_corr_length/2));
SC3_prn = [0 Bpsk_one Bpsk_zero];

[SC1_corr lag_1] = xcorr(SC1_prn, 'coeff');
[SC2_corr lag_2] = xcorr(SC2_prn, 'coeff');
[SC3_corr lag_3] = xcorr(SC3_prn, 'coeff');

index_1 = 2*lag_1/length(lag_1);
index_2 = 2*lag_2/length(lag_2);
index_3 = 2*lag_3/length(lag_2);

figure(2)
subplot(1,2,1);
plot(index_1, SC1_corr, index_3, SC3_corr, 'r--');
xlabel('��Ƭ');
ylabel('��һ������غ���');
axis([-1 1 -0.6 1]);
legend('BOC(1,1)','BPSK-R(2)');
title('BOC(1,1)�ź���BPSK�ź���������ܶԱ�');
subplot(1,2,2);
plot(index_1, SC1_corr, index_2, SC2_corr, 'r--');
xlabel('��Ƭ');
ylabel('��һ������غ���');
legend('BOC(1,1)','BOC(6,1)');
title('BOC�ź�����غ���');

figure(3)
plot(index_2, SC2_corr, index_3, SC3_corr, 'r--');
xlabel('��Ƭ');
ylabel('��һ������غ���');
legend('BOC(6,1)','BPSK-R(2)');
title('BOC�ź���������ܷ���');

%%MBOC�źŹ������ܶȷ���
S_cboc = 10/11*S1_pow + 1/11*S2_pow;

figure(4)
subplot(2,2,1);
plot(S1_pow);
xlabel('f/kHz');
ylabel('Amplitude');
title('BOC(1,1)');
subplot(2,2,2);
plot(S2_pow);
xlabel('f/kHz');
ylabel('Amplitude');
title('BOC(6,1)');
subplot(2,2,3);
plot(S_cboc);
xlabel('f/kHz');
ylabel('Amplitude');
title('MBOC��6,1,1/11��');
subplot(2,2,4);
plot(S_cboc,'r');
hold on;
plot(S1_pow);
xlabel('f/kHz');
ylabel('Amplitude');
legend('MBOC��6,1,1/11��','BOC(1,1)');

%%MBOC�ź���������ܷ���
SC_mboc_corr = 10/11*SC1_corr + 1/11*SC2_corr;

figure(5)
subplot(2,2,1);
plot(index_1, SC1_corr, index_3, SC3_corr, 'r--');
xlabel('��Ƭ');
ylabel('��һ������غ���');
axis([-1 1 -0.6 1]);
legend('BOC(1,1)','BPSK-R(2)');
title('BOC(1,1)�ź���BPSK�ź���������ܶԱ�');
subplot(2,2,2);
plot(index_2, SC2_corr, index_3, SC3_corr, 'r--');
xlabel('��Ƭ');
ylabel('��һ������غ���');
legend('BOC(6,1)','BPSK-R(2)');
title('BOC(6,1)�ź���BPSK�ź���������ܶԱ�');
subplot(2,2,3);
plot(index_1, SC1_corr, 'r--', index_2, SC2_corr);
xlabel('��Ƭ');
ylabel('��һ������غ���');
legend('BOC(1,1)','BOC(6,1)');
title('BOC�źŵ�����غ���');
subplot(2,2,4);
plot(index_1, SC1_corr,'r--');
hold on;
plot(index_1, SC_mboc_corr);
xlabel('��Ƭ');
ylabel('��һ������غ���');
axis([-1 1 -0.6 1]);
legend('BOC(1,1)','MBOC(6,1,1/11)');
title('BOC��MBOC�źŵ�����غ���');

%%GPS Galileo BDS��MBOC�ź���������ܶԱ�
SC1_SC2_corr = xcorr(SC1_prn, SC2_prn, 'coeff'); %���������

SC_Tmboc_corr = 29/33*SC1_corr + 4/33*SC2_corr;  %GPS L1C�ź�
SC_Qmboc_corr = 29/33*SC1_corr + 4/33*SC2_corr;  %BDS B1C�ź�
SC_Cboc1_corr = 10/11*SC1_corr + 1/11*SC2_corr + 2*sqrt(10/11*1/11)*SC1_SC2_corr; %Galileo E1�ź�
SC_Cboc2_corr = 10/11*SC1_corr + 1/11*SC2_corr - 2*sqrt(10/11*1/11)*SC1_SC2_corr; %Galileo E1�ź�

figure(6)
plot(index_1, SC1_corr,'k:');           %BOC(1,1)��Ϊ�ο�
hold on;
plot(index_1, SC_Tmboc_corr, 'y-');     %TMBOC�ź���QMBOC�ź�����غ���һ�������Ǵ���
hold on;
plot(index_1, SC_Qmboc_corr, 'b-');     %QMBOC
hold on;
plot(index_1, SC_Cboc1_corr, 'm--');    %CBOC+
hold on;
plot(index_1, SC_Cboc2_corr, 'r-.');    %CBOC-
xlabel('��Ƭ');
ylabel('��һ������غ���');
axis([-1 1 -0.6 1]);
legend('BOC(1,1)','TMBOC(6,1,4/33)','QMBOC(6,1,4/33)','CBOC(6,1,1/11,+)','CBOC(6,1,1/11,-)');
title('����ϵͳ�źŵ���������ܶԱ�');

%%-------------------------------------------------------------------------
%�������ڼ����Ӳ��ַ������
SB1C = 1/2*PseduCode.*SC1.*cos(2*pi*f_medium*t)+...
    sqrt(1/11)*PseduCode.*SC2.*cos(2*pi*f_medium*t)+...
    sqrt(-1)*sqrt(29/44)*PseduCode.*SC1.*cos(2*pi*f_medium*t);

SB1C_fft = fft(SB1C);
SB1C_pow = abs(SB1C_fft).^2;

figure(7)
plot(10*log10(SB1C_pow));
xlabel('Frequency(kHz)');
ylabel('PSD(dBw)');
%title('����B1C�źŵĹ������ܶ�');
grid on;
axis([0 3.5e4 -10 70]);