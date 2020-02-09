%Author:LSQ
%Date:2019/4
%Description: ����2.1.2�½�.

clc;
close all;

set(0,'defaultfigurecolor','w'); %������ͼ��������Ϊ��ɫ

%%����B1C�źŵĹ������ܶ�
load('C.mat');
load('W.mat');

f_medium = 24.58e6;      %��ƵƵ��
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

SB1C = 1/2*PseduCode.*SC1.*cos(2*pi*f_medium*t)+...
    sqrt(1/11)*PseduCode.*SC2.*cos(2*pi*f_medium*t)+...
    sqrt(-1)*sqrt(29/44)*PseduCode.*SC1.*cos(2*pi*f_medium*t);

SB1C_fft = fft(SB1C);
SB1C_pow = abs(SB1C_fft).^2;

figure(1)
plot(10*log10(SB1C_pow/f_sample));
xlabel('Ƶ��(MHz)');
ylabel('�������ܶ�(dBw/Hz)');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([0 7e4 -90 -10]);
set(gca,'xticklabel',{'0','10','20','30','40','50','60','70'});%ת��ΪMHz

%%����BOC�����źŹ������ܶ�,ֱ�Ӳ������Ϲ�ʽ��ͼ
alfa = 1;
alfa2 = 6;
beta = 1;
beta2 = 1;
beta_bpsk = 1;
pi = 3.1415926535898;
j = sqrt(-1);
fs = alfa*f0;
fs2 = alfa2*f0;
n = 2*alfa/beta;         %BOC(1,1)���ƽ���
n2 = 2*alfa2/beta2;      %BOC(6,1)���ƽ���
ts = 1/(2*fs);
ts2 = 1/(2*fs2);
tc = 1./(beta.*f0);
tc2 = 1./(beta2.*f0);
tc_bpsk = 1./(beta_bpsk.*f0);
f = -20./tc:0.001./tc:20./tc;
f_2 = -20./tc2:0.001./tc2:20./tc2;
f_bpsk = -20./tc_bpsk:0.001./tc_bpsk:20./tc_bpsk;

f_x_boc = n.*ts.*sinc(f.*ts).^2.*sin(pi.*f.*ts).^2;
f_x_bpsk = tc_bpsk.*sinc(f_bpsk.*tc_bpsk).^2;
f_x_boc2 = 1./ts2./n2.*(sin(pi*f_2.*ts2).*sin(pi*f_2.*n2.*ts2)./f_2/pi./cos(pi*f_2.*ts2)).^2;
f_mboc = 29./33.*f_x_boc + 4./33.*f_x_boc2;
f_qmboc = 29./33.*f_x_boc + j.*4./33.*f_x_boc2;

figure(2)
plot(f_bpsk,10.*log10(f_x_bpsk),'r-');axis([-8/tc_bpsk,8./tc_bpsk,-110,-60]);
hold on
%plot(f,10.*log10(f_x_boc),'k-');axis([-8./tc,8./tc,-110,-60]);
plot(f,10.*log10(f_x_boc),'b-');axis([-8./tc,8./tc,-110,-60]);
%hold on
%plot(f_2,10.*log10(f_x_boc2),'b-');axis([-8./tc2,8./tc2,-110,-60]);
hold on
plot(f_2,10.*log10(f_qmboc),'g-');axis([-8./tc2,8./tc2,-110,-60]);
xlabel('Ƶ��(MHz)');
ylabel('�������ܶ�(dBW/Hz)');
%legend('BPSK','BOC(1,1)','BOC(6,1)','QMBOC(6,1,4/33)');
legend('BPSK','BOC(1,1)','QMBOC(6,1,4/33)');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xticklabel',{'-8','-6','-4','-2','0','2','4','6','8'});%ת��ΪMHz

%%��������غ���
%BOC�ź���������ܷ���
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

%MBOC�ź���������ܷ���
SC_mboc_corr = 10/11*SC1_corr + 1/11*SC2_corr;

figure(3)
subplot(2,2,1);
plot(index_1, SC1_corr, index_3, SC3_corr, 'r--');
xlabel('��λ��ʱ(��Ƭ)');
ylabel('��һ������غ���');
axis([-1 1 -0.6 1]);
legend('BOC(1,1)','BPSK');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
title('BOC(1,1)�ź���BPSK�ź���������ܶԱ�');
subplot(2,2,2);
plot(index_2, SC2_corr, index_3, SC3_corr, 'r--');
xlabel('��λ��ʱ(��Ƭ)');
ylabel('��һ������غ���');
legend('BOC(6,1)','BPSK');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
title('BOC(6,1)�ź���BPSK�ź���������ܶԱ�');
subplot(2,2,3);
plot(index_1, SC1_corr, 'r--', index_2, SC2_corr);
xlabel('��λ��ʱ(��Ƭ)');
ylabel('��һ������غ���');
legend('BOC(1,1)','BOC(6,1)');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
title('BOC�źŵ�����غ���');
subplot(2,2,4);
plot(index_1, SC1_corr,'r--');
hold on;
plot(index_1, SC_mboc_corr);
xlabel('��λ��ʱ(��Ƭ)');
ylabel('��һ������غ���');
axis([-1 1 -0.6 1]);
legend('BOC(1,1)','MBOC(6,1,1/11)');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
title('BOC��MBOC�źŵ�����غ���');

%GPS Galileo BDS��MBOC�ź���������ܶԱ�
SC1_SC2_corr = xcorr(SC1_prn, SC2_prn, 'coeff'); %���������

SC_Tmboc_corr = 29/33*SC1_corr + 4/33*SC2_corr;  %GPS L1C�ź�
SC_Qmboc_corr = 29/33*SC1_corr + 4/33*SC2_corr;  %BDS B1C�ź�
SC_Cboc1_corr = 10/11*SC1_corr + 1/11*SC2_corr + 2*sqrt(10/11*1/11)*SC1_SC2_corr; %Galileo E1�ź�
SC_Cboc2_corr = 10/11*SC1_corr + 1/11*SC2_corr - 2*sqrt(10/11*1/11)*SC1_SC2_corr; %Galileo E1�ź�

figure(4)
plot(index_1, SC1_corr,'k:');           %BOC(1,1)��Ϊ�ο�
hold on;
plot(index_1, SC_Tmboc_corr, 'y-');     %TMBOC�ź���QMBOC�ź�����غ���һ�������Ǵ���
hold on;
plot(index_1, SC_Qmboc_corr, 'b-');     %QMBOC
hold on;
plot(index_1, SC_Cboc1_corr, 'm--');    %CBOC+
hold on;
plot(index_1, SC_Cboc2_corr, 'r-.');    %CBOC-
xlabel('��λ��ʱ(��Ƭ)');
ylabel('��һ������غ���');
axis([-1 1 -0.6 1]);
legend('BOC(1,1)','TMBOC(6,1,4/33)','QMBOC(6,1,4/33)','CBOC(6,1,1/11,+)','CBOC(6,1,1/11,-)');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);


