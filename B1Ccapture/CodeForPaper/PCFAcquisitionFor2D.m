%Author:LSQ
%Date:2019/4
%Description: Ϊѧλ��������Աȵ�PCF�����㷨�������˲�����2Dͼ.

clc;
close all;

set(0,'defaultfigurecolor','w'); %������ͼ��������Ϊ��ɫ

%%������Ա���B1C�źŵĵ�Ƶ�������в���
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
code_data_r = generatedatacode(2);  %�����źŵ����ݷ���
codeSample_r = code_r(mod(floor(t*Rc),10230)+1);
codeSample_data_r = code_data_r(mod(floor(t*Rc),10230)+1);
Qmboc_p = sqrt(1/11)*codeSample_r.*subcarr2 + ...
    j*sqrt(29/44)*codeSample_r.*subcarr1;
Boc_d = 1/2*codeSample_data_r.*subcarr1;    %Ĭ�ϴ���ʱ���ڵĵ�������ȫΪ1

Signal_B1C = Qmboc_p + Boc_d;

code_sample = floor(f_sample/Rc);   %������Ƭ����Ӧ�Ĳ�����
num_signal = length(Signal_B1C);
delay = 306*code_sample;            %��α���趨����λ��ʱ
Signal_B1C_delay = [Signal_B1C(delay : num_signal) Signal_B1C(1 : delay-1)];

IF = 24.58e6;     %��ƵƵ��
fd = 1200;        %������Ƶ��
signal_IF = Signal_B1C_delay.*cos(2*pi*(IF+fd)*t); %ģ����Ƶ�źţ�ֻ����IQ������I����

signal = awgn(signal_IF, -25);    %�Ӹ�˹������

%%����PCF�ı���B1C�źŲ����㷨
FdSearchStep = 50;      %[Hz]
DopplerRange = 5000;      %[Hz]

FdVect= -DopplerRange:FdSearchStep:DopplerRange;     %������Ƶ��������Χ

%�������ز��������
prn_p = generatecode(2);
prn_d = generatedatacode(2);
index_code = mod(floor(Rc*t),10230)+1;
prn_local = prn_p(index_code);
prn_d_local = prn_d(index_code);

%%��Ƶ�ź�QMBOC(6,1,4/33)��BOC(1,1)
idx1 = mod(floor(12*Rc*t),12)+1;
prn1_qmboc11 = [j*sqrt(20/44),j*sqrt(20/44),j*sqrt(20/44),j*sqrt(20/44),j*sqrt(20/44),...
    j*sqrt(20/44),0,0,0,0,0,0];
s1_qmboc11 = prn1_qmboc11(idx1).*prn_local;
[g1_qmboc11 x]= xcorr(Qmboc_p, s1_qmboc11, 'coeff');
prn12_qmboc11 = [0,0,0,0,0,0,j*sqrt(20/44),j*sqrt(20/44),j*sqrt(20/44),j*sqrt(20/44),...
    j*sqrt(20/44),j*sqrt(20/44)];
s12_qmboc11 = prn12_qmboc11(idx1).*prn_local;
g12_qmboc11 = xcorr(Qmboc_p, s12_qmboc11, 'coeff');

corr_sum_qmboc11 = abs(g1_qmboc11)+abs(g12_qmboc11)-abs(g1_qmboc11+g12_qmboc11);
corr_qmboc11 = xcorr(Qmboc_p, Qmboc_p, 'coeff');

%������߶ȱ任
index = x / floor(f_sample/Rc);

%��ͼ��Ƶ����
figure(1)
plot(index, corr_qmboc11,'b',index,corr_sum_qmboc11,'m');
xlabel('��Ƭ');
ylabel('��һ������غ���');
legend('ACF','PCF');
axis([-1.5 1.5 -0.5 2]);

%%����ΪPCF�����㷨��֤
%���ɾ������ڴ���ؽ��
C = zeros(length(FdVect),Num_int);     %����������Ƭ����ؽ��
C_t = zeros(1,Num_int);     %���ڴ�����ֵ
idx = 1;     %��������

for ind_FD= 1:length(FdVect)
    corr_temp = zeros(1,Num_int) ;
    fd_ind = FdVect(ind_FD);
    %�����ز�
    m = 1:Num_int;
    carrI = cos(2*pi*(IF+fd_ind)*m/f_sample);
    carrQ = sin(2*pi*(IF+fd_ind)*m/f_sample);
    for M = 0 : (Non_Coh_Sums - 1)
    %�±�Ƶ
    SigIN = signal(1+M*Num_int : Num_int+M*Num_int);
    SigOUTI = SigIN .* carrI;
    SigOUTQ = SigIN .* carrQ;
    %���ص�Ƶ��
    S1_qmboc11 = s1_qmboc11(1:Num_int);
    S12_qmboc11 = s12_qmboc11(1:Num_int);
    
    PRNLOCFFT_boc11_E = conj(fft(S1_qmboc11));
    PRNLOCFFT_boc11_L = conj(fft(S12_qmboc11));
    
    %�Ի����źŽ���FFT
    SigOUT = SigOUTI + SigOUTQ;
    Signal_fft = fft(SigOUT);
    
    %�ع���غ���
    R_boc_prn_E_11 = Signal_fft.*PRNLOCFFT_boc11_E;
    R_boc_prn_L_11 = Signal_fft.*PRNLOCFFT_boc11_L;
    
    R_E_11 = ifft(R_boc_prn_E_11);
    R_L_11 = ifft(R_boc_prn_L_11);

    R_EL_11 = R_E_11 + R_L_11;
    
    corr_temp = corr_temp + abs(R_E_11) + abs(R_L_11) - abs(R_EL_11);
    end
       
    C(idx,:) = corr_temp;
    idx = idx + 1;
end

[value, ind_mixf] = max(max(C'));
[value, ind_mixc] = max(max(C));

code_phase = (Num_int-ind_mixc)/code_sample;
doppler =(ind_mixf-1)*FdSearchStep - DopplerRange;   %[HZ]

Threshold = 3.72*10^(2.5)/0.01;    %�õ��̶�����ֵ

for i=1:Num_int
     C_t(i) = Threshold;
end

if C(ind_mixf,ind_mixc) > Threshold
    data = sprintf('The acquisition result\n Code phase:%f ��Ƭ\nDoppler frequency:%f Hz\nThreshold:%f \n',...
        code_phase,doppler,Threshold);
    disp(data);
else
    data = sprintf('Acquisition failed!\n');
    disp(data);
end

%����άͼ
[C_y, C_x]=size(C);    %C_xΪ����λ��C_yΪ������Ƶ��
X=1:C_x;Y=1:C_y;       %XΪ����λ��YΪ������Ƶ��
[x,y]=meshgrid(X,Y);
C_x_2D = (doppler+DopplerRange)/FdSearchStep+1;
C_2D = C(C_x_2D,:);
figure(2)
plot((Num_int-x)/code_sample, C_2D,'b',(Num_int-x)/code_sample, C_t,'r--');  %��������ֵ
legend('������Ƶ�ƣ�1200Hz');
axis([code_phase-2 code_phase+2 0 value+1e4]);
%axis([0 10230 0 value+1e4]);
xlabel('����λ��ʱ(��Ƭ)');ylabel('���ֵ');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

C_PCF = C;

