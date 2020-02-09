%Author:LSQ
%Date:2019/4
%Description: ���ı���B1Cģ����Ƶ�źŵ�����.

clc;
close all;

set(0,'defaultfigurecolor','w'); %������ͼ��������Ϊ��ɫ

%%���ɱ���B1Cģ����Ƶ�ź�
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

f_sample = 40*1.023e6;             %����Ƶ��
f_sc_a = 1.023e6 ;                 %���ݷ������ز�����
f_sc_b = 6*1.023e6 ;               %��Ƶ�������ز�����

t = 0 : 1/f_sample : 10e-3;

j=sqrt(-1);

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

IF = 24.58e6;     %��ƵƵ��
fd = 1500;        %������Ƶ��
Signal_B1C_out = Signal_B1C.*(cos(2*pi*(IF+fd)*t) + j*sin(2*pi*(IF+fd)*t)); %ģ����Ƶ�źţ�����IQ����

data = awgn(Signal_B1C_out, -20);    %�Ӹ�˹������

Signal_pow = abs(fft(Signal_B1C_out)).^2;
Signal_snr_pow = abs(fft(data)).^2;

figure(1)
subplot(2,1,1);
plot(t,Signal_B1C_out);
xlabel('ʱ��(ms)');
ylabel('��ֵ(v)');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([0,0.01,-1.5,1.5]);
set(gca,'xticklabel',{'0','1','2','3','4','5','6','7','8','9','10'});%ת��Ϊms
title('����B1C��Ƶ�ź�ʱ����');

subplot(2,1,2);
plot(10*log10(Signal_pow/f_sample));
xlabel('Ƶ��(MHz)');
ylabel('������(dBw/Hz)');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([0,4e5,-90,10]);
set(gca,'xticklabel',{'0','5','10','15','20','25','30','35','40'});%ת��ΪMHz
title('����B1C��Ƶ�źŹ������ܶ�');

figure(2)
subplot(2,1,1);
plot(t,data);
xlabel('ʱ��(ms)');
ylabel('��ֵ(v)');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([0,0.01,-30,30]);
set(gca,'xticklabel',{'0','1','2','3','4','5','6','7','8','9','10'});%ת��Ϊms
title('����B1C��Ƶ�ź�ʱ����');

subplot(2,1,2);
plot(10*log10(Signal_snr_pow/f_sample));
xlabel('Ƶ��(MHz)');
ylabel('������(dBw/Hz)');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([0,4e5,-60,40]);
set(gca,'xticklabel',{'0','5','10','15','20','25','30','35','40'});%ת��ΪMHz
title('����B1C��Ƶ�źŹ������ܶ�');

