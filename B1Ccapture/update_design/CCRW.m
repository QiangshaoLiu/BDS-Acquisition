clc;
close all;

%%��������Ա���B1C�źŵĵ�Ƶ�������в���
%�����������
f_sample = 36*1.023e6;             %����Ƶ��
f_sc_a = 1.023e6 ;                 %���ݷ������ز�����
f_sc_b = 6*1.023e6 ;               %��Ƶ�������ز�����
Rc = 1.023e6;                      %����������
T_process = 25e-3;                 %����ʱ��
T_int = 10e-3;                      %�������ʱ��
Non_Coh_Sums = 2;                  %(Non_Coh_Sums*T_int)ms����ɻ���ʱ��
t = 0 : 1/f_sample : T_process - 1/f_sample;
n = 0:length(t)-1;                 
j=sqrt(-1);
pi = 3.141592654;                  %Բ����
Num_int = floor(f_sample * T_int); %��ɻ���ʱ������Ӧ�Ĳ�������

%%ģ����������ź�
subcarr1 = sign(sin(2*pi*f_sc_a*t));
subcarr1(1) = 1;
subcarr2 = sign(sin(2*pi*f_sc_b*t));
subcarr2(1) = 1;
code_r = generatecode(2);           %�����ź���PRN=2����Ƶ�����е���
codeSample_r = code_r(mod(floor(t*Rc),10230)+1);
Qmboc_p = sqrt(1/11)*codeSample_r.*subcarr2 + ...
    j*sqrt(29/44)*codeSample_r.*subcarr1;

BOC_6_1 = codeSample_r.*subcarr2;
BOC_1_1 = codeSample_r.*subcarr1;

code_sample = floor(f_sample/Rc);   %������Ƭ����Ӧ�Ĳ�����
num_boc = length(Qmboc_p);
delay = 306*code_sample;            %��α���趨����λ��ʱ
Qmboc_delay = [Qmboc_p(delay : num_boc) Qmboc_p(1 : delay-1)];

IF = 24.58e6;     %��ƵƵ��
fd = 1240;        %������Ƶ��
signal_p = Qmboc_delay.*cos(2*pi*(IF+fd)*t); %ģ����Ƶ�źţ�ֻ����IQ������I����

signal = awgn(signal_p, -25);    %�Ӹ�˹������

%%����ASPeCT�ı���B1C�źŲ����㷨
FdSearchStep = 40;      %[Hz]
DopplerRange = 5000;      %[Hz]

FdVect= -DopplerRange:FdSearchStep:DopplerRange;     %������Ƶ��������Χ

%�������ز��������
prn_p = generatecode(2);
index_code = mod(floor(Rc*t),10230)+1;
prn_local = prn_p(index_code);

%%��Ƶ�ź�QMBOC(6,1,4/33)
idx1 = mod(floor(12*Rc*t),12)+1;
prn1_qmboc11 = [j*sqrt(29/44),0,0,0,0,0,0,0,0,0,0,0];
s1_qmboc11 = prn1_qmboc11(idx1).*prn_local;
[g1_qmboc11 x]= xcorr(Qmboc_p, s1_qmboc11, 'coeff');
prn1_qmboc61 = [sqrt(1/11),sqrt(1/11),sqrt(1/11),sqrt(1/11),...
    sqrt(1/11),sqrt(1/11),0,0,0,0,0,0];
s1_qmboc61 = prn1_qmboc61(idx1).*prn_local;
g1_qmboc61 = xcorr(Qmboc_p, s1_qmboc61, 'coeff');
prn12_qmboc11 = [0,0,0,0,0,0,0,0,0,0,0,j*sqrt(29/44)];
s12_qmboc11 = prn12_qmboc11(idx1).*prn_local;
g12_qmboc11 = xcorr(Qmboc_p, s12_qmboc11, 'coeff');
prn12_qmboc61 = [0,0,0,0,0,0,sqrt(1/11),sqrt(1/11),sqrt(1/11),sqrt(1/11),sqrt(1/11),sqrt(1/11)];
s12_qmboc61 = prn12_qmboc61(idx1).*prn_local;
g12_qmboc61 = xcorr(Qmboc_p, s12_qmboc61, 'coeff');

corr_sum_qmboc = abs(g1_qmboc11)+abs(g1_qmboc61)+abs(g12_qmboc11)+abs(g12_qmboc61)...
    -abs(g1_qmboc11+g12_qmboc11)-abs(g1_qmboc61+g12_qmboc61);
max_num = max(corr_sum_qmboc);
corr_qmboc = xcorr(Qmboc_p, Qmboc_p, 'coeff'); 

%%BOC(1,1)�ź�����ر߷�����
prn1_boc11 = [1,1,1,1,1,1,0,0,0,0,0,0];
s1_boc11 = prn1_boc11(idx1).*prn_local;
g1_boc11 = xcorr(BOC_1_1, s1_boc11, 'coeff');

prn2_boc11 = [0,0,0,0,0,0,1,1,1,1,1,1];
s2_boc11 = prn2_boc11(idx1).*prn_local;
g2_boc11 = xcorr(BOC_1_1, s2_boc11, 'coeff');

corr_sum_boc11 = abs(g1_boc11)+abs(g2_boc11)-abs(g1_boc11+g2_boc11);
corr_boc11 = xcorr(BOC_1_1, BOC_1_1, 'coeff'); 

%%BOC(6,1)�ź�����ر߷�����
prn1_boc61 = [1,0,0,0,0,0,0,0,0,0,0,0];
s1_boc61 = prn1_boc61(idx1).*prn_local;
g1_boc61 = xcorr(BOC_6_1, s1_boc61, 'coeff');
prn2_boc61 = [0,0,0,0,0,0,0,0,0,0,0,1];
s2_boc61 = prn2_boc61(idx1).*prn_local;
g2_boc61 = xcorr(BOC_6_1, s2_boc61, 'coeff');
prn3_boc61 = [0,1,0,0,0,0,0,0,0,0,0,0];
s3_boc61 = prn3_boc61(idx1).*prn_local;
g3_boc61 = xcorr(BOC_6_1, s3_boc61, 'coeff');
prn4_boc61 = [0,0,0,0,0,0,0,0,0,0,1,0];
s4_boc61 = prn4_boc61(idx1).*prn_local;
g4_boc61 = xcorr(BOC_6_1, s4_boc61, 'coeff');
prn5_boc61 = [0,0,1,0,0,0,0,0,0,0,0,0];
s5_boc61 = prn5_boc61(idx1).*prn_local;
g5_boc61 = xcorr(BOC_6_1, s5_boc61, 'coeff');
prn6_boc61 = [0,0,0,0,0,0,0,0,0,1,0,0];
s6_boc61 = prn6_boc61(idx1).*prn_local;
g6_boc61 = xcorr(BOC_6_1, s6_boc61, 'coeff');
prn7_boc61 = [0,0,0,1,0,0,0,0,0,0,0,0];
s7_boc61 = prn7_boc61(idx1).*prn_local;
g7_boc61 = xcorr(BOC_6_1, s7_boc61, 'coeff');
prn8_boc61 = [0,0,0,0,0,0,0,0,1,0,0,0];
s8_boc61 = prn8_boc61(idx1).*prn_local;
g8_boc61 = xcorr(BOC_6_1, s8_boc61, 'coeff');
prn9_boc61 = [0,0,0,0,1,0,0,0,0,0,0,0];
s9_boc61 = prn9_boc61(idx1).*prn_local;
g9_boc61 = xcorr(BOC_6_1, s9_boc61, 'coeff');
prn10_boc61 = [0,0,0,0,0,0,0,1,0,0,0,0];
s10_boc61 = prn10_boc61(idx1).*prn_local;
g10_boc61 = xcorr(BOC_6_1, s10_boc61, 'coeff');
prn11_boc61 = [0,0,0,0,0,1,0,0,0,0,0,0];
s11_boc61 = prn11_boc61(idx1).*prn_local;
g11_boc61 = xcorr(BOC_6_1, s11_boc61, 'coeff');
prn12_boc61 = [0,0,0,0,0,0,1,0,0,0,0,0];
s12_boc61 = prn12_boc61(idx1).*prn_local;
g12_boc61 = xcorr(BOC_6_1, s12_boc61, 'coeff');

% corr_sum_boc61 = abs(g1_boc61)+abs(g2_boc61)+abs(g3_boc61)+abs(g4_boc61)+...
%     abs(g5_boc61)+abs(g6_boc61)+abs(g7_boc61)+abs(g8_boc61)+...
%     abs(g9_boc61)+abs(g10_boc61)+abs(g11_boc61)+abs(g12_boc61)-...
%     abs(g1_boc61+g2_boc61)-abs(g3_boc61+g4_boc61)-abs(g5_boc61+g6_boc61)-abs(g7_boc61+g8_boc61)-...
%     abs(g9_boc61+g10_boc61)-abs(g11_boc61+g12_boc61);
corr_sum_boc61 = abs(g1_boc61)+abs(g2_boc61) - abs(g1_boc61+g2_boc61);
corr_boc61 = xcorr(BOC_6_1, BOC_6_1, 'coeff'); 

%������߶ȱ任
index = x / floor(f_sample/Rc);
 
 figure(1)
 plot(index, corr_qmboc,'b',index,corr_sum_qmboc,'m');
 xlabel('��Ƭ');
 ylabel('��һ������غ���');
 axis([-1.5 1.5 -0.5 2.1]);
 
 figure(2)
 plot(index, corr_sum_boc11,'b',index,corr_boc11,'m');
 xlabel('��Ƭ');
 ylabel('��һ������غ���');
 axis([-1.5 1.5 -0.6 1.5]);
 
 figure(3)
 plot(index,corr_boc61,'b',index, corr_sum_boc61,'r');
 legend('BOC(6,1)�����','�ع���غ���');
 xlabel('��Ƭ');
 ylabel('��һ������غ���');
 axis([-1.5 1.5 -1 1]);
 
 bocprncorr = xcorr(BOC_1_1,prn_local,'coeff');
 figure(4)
 plot(index, g1_boc11,'b',index,g2_boc11,'r');
 legend('BOC-d1�����','BOC-d2�����');
 xlabel('��Ƭ');
 ylabel('��һ������غ���');
 axis([-1.5 1.5 -1 1]);
 
 %%���ڶԱ�ASPeCT�㷨
 ASPeCT_11 = abs(corr_boc11).^2 - abs(bocprncorr).^2;
 x = abs(g1_boc11).^2 + abs(g2_boc11).^2 - abs(g1_boc11+g2_boc11).^2;%��ASPeCT��ȫһ��
 figure(5)
 plot(index, ASPeCT_11,'b',index,corr_sum_boc11,'r');
 legend('ASPeCT�����㷨','�Ľ���Ĳ����㷨');
 xlabel('��Ƭ');
 ylabel('��һ������غ���');
 axis([-1.5 1.5 -0.5 1.5]);
 
 qmbocprncorr = xcorr(Qmboc_p,prn_local,'coeff');
 ASPeCT_qmboc = abs(corr_qmboc).^2 - abs(qmbocprncorr).^2;
 figure(6)
 plot(index, ASPeCT_qmboc,'b',index,corr_sum_qmboc/max_num,'r');
 legend('ASPeCT�����㷨','�Ľ���Ĳ����㷨');
 xlabel('��Ƭ');
 ylabel('��һ������غ���');
 axis([-1.5 1.5 -0.5 1.1]);
 
 

 