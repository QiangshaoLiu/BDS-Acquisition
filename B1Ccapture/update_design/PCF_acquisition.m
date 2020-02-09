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

corr_sum_boc61 = abs(g1_boc61)+abs(g2_boc61)-abs(g1_boc61+g2_boc61);
corr_boc61 = xcorr(BOC_6_1, BOC_6_1, 'coeff'); 

%������߶ȱ任
index = x / floor(f_sample/Rc);
 
 figure(1)
 plot(index, corr_qmboc,'b',index,corr_sum_qmboc/max_num,'m');
 xlabel('��Ƭ');
 ylabel('��һ������غ���');
 axis([-1.5 1.5 -0.5 1.1]);
 
 figure(2)
 plot(index, corr_sum_boc11,'b',index,corr_boc11,'m');
 xlabel('��Ƭ');
 ylabel('��һ������غ���');
 axis([-1.5 1.5 -0.5 2.1]);
 
 figure(3)
 plot(index, corr_sum_boc61,'b',index,corr_boc61,'m');
 xlabel('��Ƭ');
 ylabel('��һ������غ���');
 axis([-1.5 1.5 -0.5 1.1]);

%%����ΪPCF�����㷨��ASPeCT�����㷨�ĶԱ�����֤

%�������ظ��ز�������ASPeCT�����㷨
subcarr1_loc = sign(sin(2*pi*f_sc_a*t));
subcarr2_loc = sign(sin(2*pi*f_sc_b*t));

%���ɾ������ڴ���ؽ��
C = zeros(length(FdVect),Num_int);     %����������Ƭ����ؽ��
C_aspect = zeros(length(FdVect),Num_int);  %���ڴ�ASPeCT�㷨����ؽ��
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
    %������
    S1_qmboc11 = s1_qmboc11(1:Num_int);
    S1_qmboc61 = s1_qmboc61(1:Num_int);
    S12_qmboc11 = s12_qmboc11(1:Num_int);
    S12_qmboc61 = s12_qmboc61(1:Num_int);
    
    PRNLOCFFT_boc11_E = conj(fft(S1_qmboc11));
    PRNLOCFFT_boc61_E = conj(fft(S1_qmboc61));
    PRNLOCFFT_boc11_L = conj(fft(S12_qmboc11));
    PRNLOCFFT_boc61_L = conj(fft(S12_qmboc61));   
    %�Ի����źŽ���FFT
    SigOUT = SigOUTI + SigOUTQ;
    Signal_fft = fft(SigOUT);
    %�ع���غ���
    R_boc_prn_E_11 = Signal_fft.*PRNLOCFFT_boc11_E;
    R_boc_prn_L_11 = Signal_fft.*PRNLOCFFT_boc11_L;
    R_boc_prn_L_61 = Signal_fft.*PRNLOCFFT_boc61_L;
    R_boc_prn_E_61 = Signal_fft.*PRNLOCFFT_boc61_E;
    
    R_E_11 = ifft(R_boc_prn_E_11);
    R_L_11 = ifft(R_boc_prn_L_11);
    R_E_61 = ifft(R_boc_prn_E_61);
    R_L_61 = ifft(R_boc_prn_L_61);

    R_EL_11 = R_E_11 + R_L_11;
    R_EL_61 = R_E_61 + R_L_61;
    
    corr_temp =corr_temp + abs(R_E_11) + abs(R_L_11) + abs(R_E_61) + abs(R_L_61)...
        - abs(R_EL_11) - abs(R_EL_61);
    
    %ASPeCT�����㷨
    corr_temp_aspect = zeros(1,Num_int);
    %���ظ��ز�
    Subcarr1_tot = subcarr1_loc(Num_int*M + 1 : Num_int*M + Num_int);
    Subcarr2_tot = subcarr2_loc(Num_int*M + 1 : Num_int*M + Num_int);
    prn_aspect = prn_local(1+M*Num_int : Num_int+M*Num_int);
    
    corr_temp_aspect = corr_temp_aspect + ...
        abs(ifft(Signal_fft.*conj(fft((j*sqrt(29/44)*Subcarr1_tot + ...
        sqrt(1/11)*Subcarr2_tot).*prn_aspect)))) -...
        abs(ifft(Signal_fft.*conj(fft(prn_aspect))));
    end
       
    C(idx,:) = corr_temp;
    C_aspect(idx,:) = sqrt(corr_temp_aspect.^2);
    idx = idx + 1;
end

[value, ind_mixf] = max(max(C'));
[value, ind_mixc] = max(max(C));

code_phase = (Num_int-ind_mixc)/code_sample;
doppler =(ind_mixf-1)*FdSearchStep - DopplerRange;   %[HZ]

%����Ӧ�������о�
Num_code = 12;                %�������Ƭ��Ԫ��Χ�Ĳο���Ƭ��Ԫ��Ŀ
ThresholdFactor = 9.34;      %�龯��Ϊ10e-6�����ޱ�������
Z = 0;                        %���ʹ���ֵ
for i=1:6
    Z = Z + C(ind_mixf,ind_mixc+i*code_sample)+C(ind_mixf,ind_mixc-i*code_sample);    
end
Z_aver = Z/Num_code;

Threshold = ThresholdFactor*Z_aver;    %�õ�����Ӧ����ֵ

if C(ind_mixf,ind_mixc) > Threshold
    data = sprintf('The acquisition result\n Code phase:%f ��Ƭ\nDoppler frequency:%f Hz\nThreshold:%f \n',...
        code_phase,doppler,Threshold);
    disp(data);
else
    data = sprintf('Acquisition failed!\n');
    disp(data);
end

%����άͼ
% x1 = Num_int-2046*code_sample;
% x2 = Num_int;
% y1 = (-2000+DopplerRange)/FdSearchStep+1;
% y2 = (2000+DopplerRange)/FdSearchStep+1;
% C_part=C(y1:y2,x1:x2);
% [C_y, C_x]=size(C_part);    %C_xΪ����λ��C_yΪ������Ƶ��
% X=1:C_x;Y=1:C_y;       %XΪ����λ��YΪ������Ƶ��
% [x,y]=meshgrid(X,Y);
%figure(4)
%mesh((Num_int-x)/code_sample,(y-1)*FdSearchStep - DopplerRange,C_part);  %����任
%mesh((x2-x1-x)/code_sample,(y+y1-2)*FdSearchStep - DopplerRange,C_part);  %��С��ʾ��Χ
%view(0,0);                        %��ά�ӽ�
% hold on
% plot3(code_phase,doppler,value,'k.','markersize',20);     %С�ڵ���
% text(code_phase,doppler,value,['X:',num2str(code_phase,4),char(10), ...
%     'Y:',num2str(doppler,4),char(10),'Z:',num2str(value),char(10)]);  %��ֵ
%xlabel('����λ/��');ylabel('������Ƶ��/Hz');zlabel('���ֵ');
%axis([0 2046 -2000 2000 0 value+1e4]);
%title('����B1C�źŲ�����');

%����άͼ
[C_y, C_x]=size(C);    %C_xΪ����λ��C_yΪ������Ƶ��
X=1:C_x;Y=1:C_y;       %XΪ����λ��YΪ������Ƶ��
%[x,y]=meshgrid(X,Y);
C_x_2D = (doppler+DopplerRange)/FdSearchStep+1;
C_2D = C(C_x_2D,:);
figure(5)
%plot((Num_int-x)/code_sample, C_2D,'b');
% axis([0 10230 0 value+1e4]);
axis([code_phase-2 code_phase+2 0 value+1e4]);
xlabel('����λ/��');ylabel('���ֵ');
% legend('������Ƶ�ƣ�1240Hz');

%%ASPeCT
%[value_aspect, ind_mixf_aspect] = max(max(C_aspect'));
%[value_aspect, ind_mixc_aspect] = max(max(C_aspect));
%code_phase_aspect = (Num_int-ind_mixc_aspect)/code_sample;
%doppler_aspect =(ind_mixf_aspect-1)*FdSearchStep - DopplerRange;   %[HZ]

%��ASPeCT�㷨����άͼ
%x1 = Num_int-2046*code_sample;
%x2 = Num_int;
%y1 = (-2000+DopplerRange)/FdSearchStep+1;
%y2 = (2000+DopplerRange)/FdSearchStep+1;
%C_part=C_aspect(y1:y2,x1:x2);
%[C_y_aspect, C_x_aspect]=size(C_part);    %C_yΪ����λ��C_xΪ������Ƶ��
%X_aspect=1:C_x_aspect;Y_aspect=1:C_y_aspect;       %XΪ����λ��YΪ������Ƶ��
%[x_aspect,y_aspect]=meshgrid(X_aspect,Y_aspect);
%figure(6)
%mesh((Num_int-x_aspect)/code_sample,(y_aspect-1)*FdSearchStep - DopplerRange,C_aspect); 
%mesh((x2-x1-x_aspect)/code_sample,(y_aspect+y1-2)*FdSearchStep - DopplerRange,C_part);  %����任
%hold on
%view(0,90);                        %��ά�ӽ�
% plot3(code_phase_aspect,doppler_aspect,value_aspect,'k.','markersize',20);     %С�ڵ���
% text(code_phase_aspect,doppler_aspect,value_aspect,['X:',num2str(code_phase_aspect,4),char(10), ...
%     'Y:',num2str(doppler_aspect,4),char(10),'Z:',num2str(value_aspect),char(10)]);  %��ֵ
%xlabel('����λ/��');ylabel('������Ƶ��/Hz');zlabel('���ֵ');
%axis([0 10230 -5000 5000 0 value_aspect+1e4]);
% axis([0 2046 -2000 2000 0 value_aspect+1e4]);

%��ASPeCT�㷨�Ķ�άͼ
% [C_y_aspect, C_x_aspect]=size(C_aspect);    %C_yΪ����λ��C_xΪ������Ƶ��
% X_aspect=1:C_x_aspect;Y_aspect=1:C_y_aspect;       %XΪ����λ��YΪ������Ƶ��
% [x_aspect,y_aspect]=meshgrid(X_aspect,Y_aspect);
% C_x_2D_aspect = (doppler_aspect+DopplerRange)/FdSearchStep+1;
% C_2D_aspect = C_aspect(C_x_2D_aspect,:);
% figure(7)
% plot((Num_int-x_aspect)/code_sample, C_2D_aspect,'b');
% axis([0 10230 0 value_aspect+1e4]);
% %%axis([code_phase_aspect-6 code_phase_aspect+6 0 value_aspect+1e4]);
% xlabel('����λ/��');ylabel('���ֵ');
% legend('������Ƶ�ƣ�1240Hz');
