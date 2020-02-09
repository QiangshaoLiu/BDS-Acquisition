%����ر߷����������ķ�������������ڼ��
%����ɻ�������ɻ������������ɻ��ֵĶԱȣ��������ڼ��

clc;
close all;

%����Weil�룬PRN��Ϊ1��10,���ݷ����͵�Ƶ��������
PRN = 1;

switch PRN
    case (1)
        phase_diff_d = 2678;
        intercept_pot_d = 699;
        phase_diff_p = 796;
        intercept_pot_p = 7575;
    case (2)
        phase_diff_d = 4802;
        intercept_pot_d = 694;
        phase_diff_p = 156;
        intercept_pot_p = 2369;
    case (3)
        phase_diff_d = 958;
        intercept_pot_d = 7318;
        phase_diff_p = 4198;
        intercept_pot_p = 5688;
    case (4)
        phase_diff_d = 859;
        intercept_pot_d = 2127;
        phase_diff_p = 3941;
        intercept_pot_p = 539;
    case (5)
        phase_diff_d = 3843;
        intercept_pot_d = 715;
        phase_diff_p = 1374;
        intercept_pot_p = 2270;
    case (6)
        phase_diff_d = 2232;
        intercept_pot_d = 6682;
        phase_diff_p = 1338;
        intercept_pot_p = 7306;
    case (7)
        phase_diff_d = 124;
        intercept_pot_d = 7850;
        phase_diff_p = 1833;
        intercept_pot_p = 6457;
    case (8)
        phase_diff_d = 4352;
        intercept_pot_d = 5495;
        phase_diff_p = 2521;
        intercept_pot_p = 6254;
    case (9)
        phase_diff_d = 1816;
        intercept_pot_d = 1162;
        phase_diff_p = 3175;
        intercept_pot_p = 5644;
    case (10)
        phase_diff_d = 1126;
        intercept_pot_d = 7682;
        phase_diff_p = 168;
        intercept_pot_p = 7119;
     otherwise
        phase_diff_d = 2678;
        intercept_pot_d = 699;
        phase_diff_p = 796;
        intercept_pot_p = 7575;    %������Χ��Ĭ��PRN=1
end 

legendrelength = 10243;
WeilCodelength = 10230;
legendre = zeros(1, legendrelength);
WeilCode_d = zeros(1, WeilCodelength);
WeilCode_p = zeros(1, WeilCodelength);

for k = 1 : legendrelength-1 
    for x = 1 : (legendrelength-1)/2
        if mod(k,legendrelength) == mod(x^2, legendrelength)
            legendre(k) = 1;
        end
    end
end

legendre = [0 legendre(1:legendrelength-1)];   %����legendre���� 

%����B1C����/��Ƶ��������
for k = 0 : WeilCodelength-1     
    WeilCode_d(k+1) = mod(sum([legendre(mod((k+intercept_pot_d-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff_d+intercept_pot_d-1), legendrelength)+1)]),2);
    WeilCode_p(k+1) = mod(sum([legendre(mod((k+intercept_pot_p-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff_p+intercept_pot_p-1), legendrelength)+1)]),2);
end

%��Weil��任�ɼ����룬0��ʾ�ߵ�ƽ��+1����1��ʾ�͵�ƽ��-1��
for k = 1 : WeilCodelength
    if WeilCode_d(k) == 0
       WeilCode_d(k) = 1;
    else
       WeilCode_d(k) = -1;
    end
    if WeilCode_p(k) == 0
       WeilCode_p(k) = 1;
    else
       WeilCode_p(k) = -1;
    end    
end

f_sample = 10*1.023e6;             %����Ƶ��
f_sc_a = 1.023e6 ;                 %���ݷ������ز�����
f_sc_b = 6*1.023e6 ;               %��Ƶ�������ز�����
T_process = 50e-3;                 %����ʱ��
T_int = 10e-3;                      %�������ʱ��
Non_Coh_Sums = 5;                  %(Non_Coh_Sums*T_int)ms����ɻ���ʱ��
t = 0 : 1/f_sample : T_process - 1/f_sample;

j=sqrt(-1);

%N = length(t);     %��������ɻ��֣�ֱ�ӶԴ���ʱ���ڵ����ݽ���ASPeCT����
N = floor(f_sample * T_int);

Subcarr1 = sign(sin(2*pi*f_sc_a*t));
Subcarr2 = sign(sin(2*pi*f_sc_b*t));

Rc = 1.023e6;       %����������
CodeSample_d = WeilCode_d(mod(floor(t*Rc),10230)+1);
CodeSample_p = WeilCode_p(mod(floor(t*Rc),10230)+1);
B1C_data = 1/2*CodeSample_d.*Subcarr1;      
B1C_poilt = sqrt(1/11)*CodeSample_p.*Subcarr2 + ...
    j*sqrt(29/44)*CodeSample_p.*Subcarr1;

%����ͨ���͵�Ƶͨ��������ź�
Signal_B1C = B1C_data + B1C_poilt;

num_sample = floor(f_sample/Rc);
num = length(Signal_B1C);
delay = 272*num_sample;       %���ź��趨��ʱ
Signal_B1C_delay = [Signal_B1C(delay : num) Signal_B1C(1 : delay-1)];

IF = 24.58e6;     %��ƵƵ��
fd = 1500;        %������Ƶ��
Signal_B1C_out = Signal_B1C_delay.*cos(2*pi*(IF+fd)*t); %ģ����Ƶ�źţ�ֻ����IQ������I����

data = awgn(Signal_B1C_out, -30);    %�Ӹ�˹������

acqSearchStep = 250;      %[Hz]
DopplerRange = 5000;      %[Hz]

FD_vect= -DopplerRange:acqSearchStep:DopplerRange;

% FFT acquisition
idx = 1;
n=0:length(t)-1;
ind_cod = mod(floor(n*Rc/f_sample),WeilCodelength)+1;
SigLOC_tot_d = WeilCode_d(ind_cod);
SigLOC_tot_p = WeilCode_p(ind_cod);

%AsPect�������õ�������غ���
%���ݷ���
[BOCPRN_corr lag] = xcorr(Subcarr1.*SigLOC_tot_d, SigLOC_tot_d, 'coeff');  %BOC�ź���PRN��Ļ���غ���
BOC_corr = xcorr(Subcarr1.*SigLOC_tot_d, 'coeff');    %PRN�������غ���
Aspect = abs(BOC_corr).^2 - abs(BOCPRN_corr).^2;
Aspect_b = abs(BOC_corr).^2 - 1.2*abs(BOCPRN_corr).^2;
%��Ƶ����
[BOCPRN_corr_p lag_p] = xcorr(Subcarr2.*SigLOC_tot_p+j*Subcarr1.*SigLOC_tot_p, SigLOC_tot_p, 'coeff');  %BOC�ź���PRN��Ļ���غ���
BOC_corr_p = xcorr(Subcarr2.*SigLOC_tot_p+j*Subcarr1.*SigLOC_tot_p, 'coeff');    %PRN�������غ���
Aspect_p = abs(BOC_corr_p).^2 - abs(BOCPRN_corr_p).^2;
Aspect_p_b = abs(BOC_corr_p).^2 - 1.2*abs(BOCPRN_corr_p).^2;

%������߶ȱ任
index = lag / num_sample;
index_p = lag_p / num_sample;

figure(1)
subplot(1,2,1);
plot(index, Aspect, 'b', index, BOC_corr, 'r--',index,Aspect_b,'g-.');
xlabel('��Ƭ');
ylabel('��һ������غ���');
axis([-1.5 1.5 -0.6 1.2]);
legend('�����߷�','ֱ������');
title('���ݷ���������ر߷�����ǰ��Ա�');
subplot(1,2,2);
plot(index_p, Aspect_p, 'b', index_p, BOC_corr_p, 'r--',index,Aspect_p_b,'g-.');
xlabel('��Ƭ');
ylabel('��һ������غ���');
axis([-1.5 1.5 -0.6 1.2]);
legend('�����߷�','ֱ������');
title('��Ƶ����������ر߷�����ǰ��Ա�');

C = zeros(length(FD_vect),N);      %��ɻ��������ɻ�������
C_non = zeros(length(FD_vect),N);  %ֻ����ɻ���

for ind_FD= 1:length(FD_vect)
    FD = FD_vect(ind_FD);
    corr = zeros(1,N) + j*zeros(1,N);
    corr_non = zeros(1,N) + j*zeros(1,N);
    %�����ز�
    m = 0:N-1;
    argx = 2*pi*(IF+FD)/f_sample;
    carrI = cos(argx*m);
    carrQ = sin(argx*m);

    %������ۼ�
    for M = 0 : (Non_Coh_Sums - 1);
    %for M = 0 : (1 - 1);
    %�±�Ƶ
    SigIN = data(N*M + 1 : N*M + N);
    SigOUTI = SigIN .* carrI;
    SigOUTQ = SigIN .* carrQ;
    %������ͱ����ز�
    PRNLOC_d = SigLOC_tot_d(N*M + 1 : N*M + N);          
    PRNLOC_p = SigLOC_tot_p(N*M + 1 : N*M + N);
    Subcarr1_tot = Subcarr1(N*M + 1 : N*M + N);
    Subcarr2_tot = Subcarr2(N*M + 1 : N*M + N);

    PRNLOC = PRNLOC_d + j*PRNLOC_p;
    PRNLOCFFT = fft(PRNLOC, N);
    PRNLOCFFT_conj = conj(PRNLOCFFT);

    BOCLOC = (Subcarr1_tot + j*Subcarr2_tot).*PRNLOC;
    BOCLOC_p = (Subcarr1_tot + j*Subcarr2_tot).*PRNLOC_p;
    BOCLOCFFT = fft(BOCLOC, N);
    BOCLOCFFT_conj = conj(BOCLOCFFT);    
    
    %ASPeCT: BOC-BOC �� BOC-PRN
    corr = corr + abs(ifft(fft(SigOUTI,N).*BOCLOCFFT_conj)).^2 - ...
        abs(ifft(fft(SigOUTI,N).*PRNLOCFFT_conj)).^2 + ...
        abs(ifft(fft(SigOUTQ,N).*BOCLOCFFT_conj)).^2 - ...
        abs(ifft(fft(SigOUTQ,N).*PRNLOCFFT_conj)).^2;
    %������ɻ���
    corr_non = abs(ifft(fft(SigOUTI,N).*BOCLOCFFT_conj)).^2 - ...
        abs(ifft(fft(SigOUTI,N).*PRNLOCFFT_conj)).^2 + ...
        abs(ifft(fft(SigOUTQ,N).*BOCLOCFFT_conj)).^2 - ...
        abs(ifft(fft(SigOUTQ,N).*PRNLOCFFT_conj)).^2;   
    end
    
    C(idx,:) = corr;
    C_non(idx,:) = corr_non;
    idx = idx+1;
end

%Find the main peak in the correlation floor and the corresponding frequency bin index
[bb, ind_mixf] = max(max(C'));
[bb, ind_mixc] = max(max(C));
[aa, ind_mixf_a] = max(max(C_non'));
[aa, ind_mixc_a] = max(max(C_non));

%code_phase = floor((int_n - ind_mixc)/num_sample);
code_phase = ceil((N - ind_mixc)/num_sample);
doppler =(ind_mixf-1)*acqSearchStep - DopplerRange;   %[HZ]
code_phase_non = ceil((N - ind_mixc_a)/num_sample);
doppler_non = (ind_mixf_a-1)*acqSearchStep - DopplerRange;

%����άͼ
[C_y, C_x]=size(C);    %C_yΪ����λ��C_xΪ������Ƶ��
X=1:C_x;Y=1:C_y;       %XΪ����λ��YΪ������Ƶ��
[x,y]=meshgrid(X,Y);
figure(2)
mesh((N-x)/num_sample,(y-1)*acqSearchStep - DopplerRange,C);  %����任
hold on
plot3(code_phase,doppler,bb,'k.','markersize',20);     %С�ڵ���
text(code_phase,doppler,bb,['X:',num2str(code_phase),char(10), ...
    'Y:',num2str(doppler),char(10),'Z:',num2str(bb),char(10)]);  %��ֵ
xlabel('����λ/��');ylabel('������Ƶ��/Hz');zlabel('���ֵ');
title('�Ľ����㷨�Ĳ�����');

[C_y_non, C_x_non]=size(C_non);    %C_yΪ����λ��C_xΪ������Ƶ��
X=1:C_x_non;Y=1:C_y_non;       %XΪ����λ��YΪ������Ƶ��
[x_non,y_non]=meshgrid(X,Y);
figure(3)
mesh((N-x_non)/num_sample,(y_non-1)*acqSearchStep - DopplerRange,C_non);  %����任
hold on
plot3(code_phase_non,doppler_non,aa,'k.','markersize',20);     %С�ڵ���
text(code_phase_non,doppler_non,aa,['X:',num2str(code_phase_non),char(10), ...
    'Y:',num2str(doppler_non),char(10),'Z:',num2str(aa),char(10)]);  %��ֵ
xlabel('����λ/��');ylabel('������Ƶ��/Hz');zlabel('���ֵ');
title('�Ľ�ǰ�㷨�Ĳ�����');

%����άͼ
C_x_2D = (doppler+DopplerRange)/acqSearchStep+1;
C_2D = C(C_x_2D,:);
C_2D_non = C_non(C_x_2D,:);

figure(4)
fuc1 = plot((N-x)/num_sample, C_2D, 'b');
hold on;
fuc2 = plot((N-x)/num_sample, C_2D_non, 'r');

legend([fuc1(1),fuc2(1)],'����ɻ�������ɻ�������','��ɻ���');
axis([code_phase-6 code_phase+6 -3e8 bb]);
xlabel('����λ/��');ylabel('���ֵ');
title('���ַ�ʽЧ���Ա�');


