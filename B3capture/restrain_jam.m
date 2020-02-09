%*********************************************************************
%������Ϊ����
%���ߣ�LSQ
%���ڣ�2016��11��24��
%---------------------------����˵��-----------------------------------
%��������б����������խ�����ŷ�����֤
%���������޷����������㷨����������ţ�Kֵ����һ�׾ط���
%��ֵ���޷���Ƶ������Ӧ���޷�����խ������
%---------------------------�������-----------------------------------
%��Ƶfm=46.52e6��������fs=62e6�������-20db�����ű�50db
%**********************************************************************
fc=1268.52e6;                        %����B3��Ƶ
fm=46.52e6;                          %��ƵƵ��
fs=62e6;                             %��Ƶ��������
Rcode=10.23e6;                       %����������
SNR=-20;                             %�����
snr=10.^(SNR/10);                    %������ϵ
ISR=50;                              %��Ƶ���ŵĸ��ű�
isr=10.^(ISR/10);                    %������ϵ
L=62000;                             %�����������Ӧ1ms����
t=0:L-1;                             %����ʱ���
sim_t=1;                             %����ʱ��1ms

%CA�����
load PNcode.mat;

%�����ź�
signal=PNcode_s.*cos(2*pi*fm*t/fs);   %��ֵΪ1
Power_sig=(signal*signal')/L;         %�źŹ���

%��˹����
Power_Guass=Power_sig/snr;            %��������������������
Guass_noise=sqrt(Power_Guass)*randn(1,L);%��˹����

%�Ե����˲������޴���20MHz,���޴����ź�ͨ����Χ��
w1=(15.48e6)-10e6;
w2=(15.48e6)+10e6;
wp=[w1*2/fs-0.01 w2*2/fs+0.01];            %20MHz���� �½�3dB����
b=fir1(100,wp);                            %100�׵Ĵ�ͨ�˲���
Guass_noise_NB=filtfilt(b,1,Guass_noise);  %�˲�
Power_noi=Guass_noise_NB*Guass_noise_NB'/L;%��������
fft_Gua_noi=fftshift(fft(Guass_noise));
fft_Gua_noi_NB=fftshift(fft(Guass_noise_NB));

%��ͼ�Աȵ�����˲��������
f=(0:L-1)*fs/L-fs/2;%Ƶ��
figure(1);
subplot(211);
plot(f,fft_Gua_noi,'b');
xlabel('f');
ylabel('fft_Gua_noi');
title('��˹����');
subplot(212);
plot(f,fft_Gua_noi_NB,'r');
xlabel('f');
ylabel('fft_Gua_noi_NB')
title('�ź�ͨ���ڵĴ��޸�˹����');

%�����źżӰ���
xif_sn=signal+Guass_noise_NB;          %�����źżӸ�˹����
figure(2);
subplot(211);
plot(t,xif_sn);                        %���յ��Ĳ������ŵ��ź�
xlabel('t');
ylabel('xif');

%����
%������ţ�ռ�ձ�Ϊ4/100�����ĸ�ʱ���ֱ�ʵʩ���ţ������ֵΪ�ź�����ֵ��400��
l_pulse=L/100;                     %ÿ�����ŵ�������ų��ȣ�620
A_pulse=400;
t0=2000:2000+l_pulse-1;
t1=12000:12000+l_pulse-1;
t2=32000:32000+l_pulse-1;
t3=52000:52000+l_pulse-1;          %�ĸ�����ʱ���
pulse=zeros(1,L);
pulse_l=randn(1,l_pulse);
pulse_l(pulse_l>0.5)=A_pulse;
pulse_l(pulse_l<0.5)=-A_pulse;
pulse(t0)=pulse_l;
pulse(t1)=pulse_l;
pulse(t2)=pulse_l;
pulse(t3)=pulse_l;                  %�����ĸ�ʱ�����������
xif_pulse=xif_sn+pulse;             %�źż��������
subplot(212);
plot(t,xif_pulse,'r');
xlabel('t');
ylabel('xif_pulse');

%ʱ������������ƣ�ʱ�������޷���ʱ���������㣬����Ҫ���޷������޵�ȷ��
%���޷�����Ϊ20������ֵ�������޵������ֱ����Ϊ16����С�ڵĲ���
%���㷨
for i=1:L
    if abs(xif_pulse(i))>20
        xif_pulse(i)=0;
    end
end
figure(3);
subplot(211);
plot(t,xif_pulse,'r');
axis([-10 L+10 -40 40]);
xlabel('t');
ylabel('xif_pulse');
title('�������㷨�����������');

%�޷���
xif_pulse=xif_sn+pulse;%�źż��������
for i=1:L
    if abs(xif_pulse(i))>20
        xif_pulse(i)=16;
    end
end
subplot(212);
plot(t,xif_pulse,'r');
axis([-10 L+10 -40 40]);
xlabel('t');
ylabel('xif_pulse');
title('�����޷��������������');

%�������źͶ������ţ���������ģ��խ������
Jam_f=[0e6,1e6,2e6,3e6];                   %��Ƶ���ž���Ƶ�ľ���
Power_jam=Power_sig*isr;
A_jam=sqrt(2*Power_jam);
Jam_1=A_jam*cos(2*pi*(fm+Jam_f(1))*t/fs);
Jam_2=A_jam*cos(2*pi*(fm+Jam_f(2))*t/fs);
Jam_3=A_jam*cos(2*pi*(fm+Jam_f(3))*t/fs);
Jam_4=A_jam*cos(2*pi*(fm+Jam_f(4))*t/fs);  %�ĸ������ź�

en=3;
if en==1                                   %enΪ����ʹ��
    Jam=Jam_1;                             %��Ƶ����
elseif en==2
    Jam=Jam_1+Jam_2;
elseif en==3;
    Jam=Jam_1+Jam_2+Jam_3;
elseif en==4;
    Jam=Jam_1+Jam_2+Jam_3+Jam_4;           %�������ţ�ģ��խ������
elseif en==0;
    Jam=zeros(1,L);                        %ʹ��Ϊ0���޸���
end 
xif_nsj=xif_sn+Jam;                        %���յ��Ĵ����ź��������ź�

%Ƶ����������㷨Ҫ��������ź���fft�任��Ϊ�˼���Ƶ��й¶�����ٸ����ź��԰�Ե����źŲ����ĸ��ţ���Ҫfft֮ǰ�ô���������ƽ������
%��������������hanning��������hamming������������blackman��������kaiser��ͨ��wȡֵѡȡ���õĴ�����
w=4;
if w==0
    window=ones(1,L);                     %���Ӵ�
elseif w==1
    window=hanning(1,L);
elseif w==2
    window=hamming(1,L);
elseif w==3
    window=blackman(1,L);
elseif w==4
    window=blackman(1,L);                 %L�㲼��������
elseif w==5
    window=kaiser(1,L);
end

xif_nsj_win=xif_nsj.*window;              %�Խ����źŴ�����ƽ������

%���������㷨���̶����ޡ���Kֵ��
%Kֵ���Ǹ��������㷨�й̶����޷���һ�֣��Ӹղŵķ�Ƶͼ�п��Կ�����û�ڸ�˹�����µĵ����źŷ�ֵ����Լ��80db���£�����ѡ���̶�����80db��Ϊ������ޣ�����80db��Ƶ�㴦ֱ�����㴦��
xif_K_fft=fft(xif_nsj_win);
f=(0:L-1)*L/fs;
xif_fft_abs=20*log10(abs(xif_K_fft));
figure(4);
subplot(211);
plot(f,xif_fft_abs);                      %���Ʒ�Ƶͼ
xlabel('f/MHz');
ylabel('xif_fft_abs');
Th_K=80;
Jam_index=find(xif_fft_abs>Th_K);         %�ҳ������źŷ�ֵ����80dB��Ƶ��
xif_K_fft(Jam_index)=0;                   %����
xif_fft_K_abs=20*log10(abs(xif_K_fft));
subplot(212);
plot(f,xif_fft_K_abs);                    %���Ʒ�Ƶͼ
xlabel('f/MHz');
ylabel('xif_fft_K_abs');
title('Kֵ������������խ������');
xif_K_ifft=ifft(xif_K_fft);               %Kֵ�������Ŵ����ifft���

%���������㷨���̶����ޡ���һ�׾ط�
%һ�׾ط�������ΪTh=q*u������uΪFFT���ֵ��qΪ�Ż�ϵ��������1024���FFT��ͨ���Ӵ��������Ƶ��й¶
xif_one_fft=fft(xif_nsj_win);%L��fft
xif_one_fft_abs=abs(xif_one_fft);
xif_one_fft_absmod=20*log10(xif_one_fft_abs);
figure(5);
subplot(211);
plot(f,xif_one_fft_absmod);               %���Ʒ�Ƶͼ
xlabel('f/MHz');
ylabel('xif_one_fft_abs');
u=mean(xif_one_fft_abs,2);                %��L����źŷ�ֵ��������mean����Ϊ����������
q=1;                                      %�Ż�ϵ������̽��ȡֵ
Th_one=q*u;                               %�������
fft_one_index=find(xif_one_fft_abs>Th_one);%�ҳ��������޵�Ƶ�㼴Ϊ��������λ��
xif_one_fft(fft_one_index)=0;              %����
xif_one_fft_abs=20*log10(abs(xif_one_fft));
subplot(212);
plot(f,xif_one_fft_abs);                   %���Ʒ�Ƶͼ
xlabel('f/MHz');
ylabel('xif_one_fft_abs');
title('һ�׾ط�����������խ������');
xif_one_ifft=ifft(xif_one_fft);            %һ�׾ؿ����Ŵ����ifft���

%���������㷨���̶����޷�����ֵ���޷�
%��ֵ���޷�Ҳ�ǹ̶����޷���һ�֣�����ȡֵΪ�����źž���fft����м�ֵ����Th=r*median{X[k]},ע�ⲻ��ƽ��ֵ�������м�ֵ��rΪ�Ż�ϵ��
xif_median_fft=fft(xif_nsj_win);           %L��fft
xif_median_fft_abs=abs(xif_median_fft);
xif_median_fft_absmod=20*log10(xif_median_fft_abs);
figure(6);
subplot(211);
plot(f,xif_median_fft_absmod);             %���Ʒ�Ƶͼ
xlabel('f/MHz');
ylabel('xif_median_fft_abs');
xif_median_fftsort=sort(xif_median_fft_abs);
median=0.5*(xif_median_fftsort(L/2)+xif_median_fftsort(L/2+1));%��fft�任������ݴ�С���������ҳ��м�ֵ
r=4;                                       %�Ż�ϵ����̽����Ϊ1.5
Th_median=r*median;                        %�������
fft_median_index=find(xif_median_fft_abs>Th_median);%�ҳ��������޵�Ƶ�㼴Ϊ��������λ��
xif_median_fft(fft_median_index)=0;        %����
xif_median_fft_abs=20*log10(abs(xif_median_fft));
subplot(212);
plot(f,xif_median_fft_abs);                %���Ʒ�Ƶͼ
xlabel('f/MHz');
ylabel('xif_median_fft_abs');
title('��ֵ���޷�����������խ������');
xif_median_ifft=ifft(xif_median_fft);      %��ֵ���޷������Ŵ����ifft���

%Ƶ���������޷�
%���յ��ĵ����źŽ���Ϊխ����˹�źţ��������������ֲ��������ƽ������ָ���ֲ�
xif_adjust_fft=fft(xif_nsj_win);%L��fft
xif_adjust_fft_abs=abs(xif_adjust_fft);
xif_adjust_fft_absmod=20*log10(xif_adjust_fft_abs);
figure(7);
subplot(211);
plot(f,xif_adjust_fft_absmod);              %���Ʒ�Ƶͼ
xlabel('f/MHz');
ylabel('xif_adjust_fft_abs');
xif_adjust_fft_abs2=abs(xif_adjust_fft).*abs(xif_adjust_fft);
u_adjust=mean(xif_adjust_fft_abs2);        %����Ƶ�ָ���ֲ��ľ�ֵ
Th_adjust=4*u_adjust;                      %����Ӧ����
%���ݰ���ƽ����ȡֵ��Χȷ��fft��ȡֵ
fft_adjust_index1=find(0*Th_adjust<xif_adjust_fft_abs2<Th_adjust);
fft_adjust_index2=find(Th_adjust<xif_adjust_fft_abs2<8*Th_adjust);
fft_adjust_index3=find(8*Th_adjust<xif_adjust_fft_abs2<64*Th_adjust);
fft_adjust_index4=find(64*Th_adjust<xif_adjust_fft_abs2<512*Th_adjust);
fft_adjust_index5=find(xif_adjust_fft_abs2>512*Th_adjust);
xif_adjust_fft(fft_adjust_index1)=xif_adjust_fft(fft_adjust_index1);
xif_adjust_fft(fft_adjust_index2)=xif_adjust_fft(fft_adjust_index2)/8;
xif_adjust_fft(fft_adjust_index3)=xif_adjust_fft(fft_adjust_index3)/64;
xif_adjust_fft(fft_adjust_index4)=xif_adjust_fft(fft_adjust_index4)/512;
xif_adjust_fft(fft_adjust_index5)=0;       %����ֵ��Ϊ0
xif_adjust_fft_abs=20*log10(abs(xif_adjust_fft));
subplot(212);
plot(f,xif_adjust_fft_abs);                %���Ʒ�Ƶͼ
xlabel('f/MHz');
ylabel('xif_adjust_fft_abs');
title('����Ӧ���޷�����������խ������');
xif_adjust_ifft=ifft(xif_adjust_fft);      %����Ӧ���޷������Ŵ����ifft���

