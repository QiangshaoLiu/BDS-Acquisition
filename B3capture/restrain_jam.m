%*********************************************************************
%化腐朽为神奇
%作者：LSQ
%日期：2016年11月24日
%---------------------------程序说明-----------------------------------
%本程序进行北斗抗脉冲和窄带干扰仿真验证
%包括脉冲限幅和脉冲置零法抑制脉冲干扰，K值法，一阶矩法，
%中值门限法，频域自适应门限法抑制窄带干扰
%---------------------------程序参数-----------------------------------
%中频fm=46.52e6，采样率fs=62e6，信噪比-20db，干信比50db
%**********************************************************************
fc=1268.52e6;                        %北斗B3载频
fm=46.52e6;                          %中频频率
fs=62e6;                             %中频采样速率
Rcode=10.23e6;                       %北斗码速率
SNR=-20;                             %信噪比
snr=10.^(SNR/10);                    %倍数关系
ISR=50;                              %单频干扰的干信比
isr=10.^(ISR/10);                    %倍数关系
L=62000;                             %仿真点数，对应1ms数据
t=0:L-1;                             %仿真时间点
sim_t=1;                             %仿真时间1ms

%CA码产生
load PNcode.mat;

%卫星信号
signal=PNcode_s.*cos(2*pi*fm*t/fs);   %幅值为1
Power_sig=(signal*signal')/L;         %信号功率

%高斯底噪
Power_Guass=Power_sig/snr;            %根据信噪比算出噪声功率
Guass_noise=sqrt(Power_Guass)*randn(1,L);%高斯底噪

%对底噪滤波将其限带在20MHz,即限带在信号通带范围中
w1=(15.48e6)-10e6;
w2=(15.48e6)+10e6;
wp=[w1*2/fs-0.01 w2*2/fs+0.01];            %20MHz带宽 下降3dB带外
b=fir1(100,wp);                            %100阶的带通滤波器
Guass_noise_NB=filtfilt(b,1,Guass_noise);  %滤波
Power_noi=Guass_noise_NB*Guass_noise_NB'/L;%噪声功率
fft_Gua_noi=fftshift(fft(Guass_noise));
fft_Gua_noi_NB=fftshift(fft(Guass_noise_NB));

%绘图对比低噪和滤波后的噪声
f=(0:L-1)*fs/L-fs/2;%频点
figure(1);
subplot(211);
plot(f,fft_Gua_noi,'b');
xlabel('f');
ylabel('fft_Gua_noi');
title('高斯低噪');
subplot(212);
plot(f,fft_Gua_noi_NB,'r');
xlabel('f');
ylabel('fft_Gua_noi_NB')
title('信号通带内的带限高斯白噪');

%导航信号加白噪
xif_sn=signal+Guass_noise_NB;          %导航信号加高斯白噪
figure(2);
subplot(211);
plot(t,xif_sn);                        %接收到的不含干扰的信号
xlabel('t');
ylabel('xif');

%干扰
%脉冲干扰，占空比为4/100，分四个时间点分别实施干扰，脉冲幅值为信号最大幅值的400倍
l_pulse=L/100;                     %每个干扰点脉冲干扰长度，620
A_pulse=400;
t0=2000:2000+l_pulse-1;
t1=12000:12000+l_pulse-1;
t2=32000:32000+l_pulse-1;
t3=52000:52000+l_pulse-1;          %四个干扰时间点
pulse=zeros(1,L);
pulse_l=randn(1,l_pulse);
pulse_l(pulse_l>0.5)=A_pulse;
pulse_l(pulse_l<0.5)=-A_pulse;
pulse(t0)=pulse_l;
pulse(t1)=pulse_l;
pulse(t2)=pulse_l;
pulse(t3)=pulse_l;                  %产生四个时间点的脉冲干扰
xif_pulse=xif_sn+pulse;             %信号加脉冲干扰
subplot(212);
plot(t,xif_pulse,'r');
xlabel('t');
ylabel('xif_pulse');

%时域脉冲干扰抑制，时域脉冲限幅和时域脉冲置零，最主要是限幅器门限的确定
%设限幅门限为20，绝对值高于门限的置零或直接置为16处理，小于的不变
%置零法
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
title('脉冲置零法抑制脉冲干扰');

%限幅法
xif_pulse=xif_sn+pulse;%信号加脉冲干扰
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
title('脉冲限幅法抑制脉冲干扰');

%单音干扰和多音干扰，多音干扰模拟窄带干扰
Jam_f=[0e6,1e6,2e6,3e6];                   %单频干扰距中频的距离
Power_jam=Power_sig*isr;
A_jam=sqrt(2*Power_jam);
Jam_1=A_jam*cos(2*pi*(fm+Jam_f(1))*t/fs);
Jam_2=A_jam*cos(2*pi*(fm+Jam_f(2))*t/fs);
Jam_3=A_jam*cos(2*pi*(fm+Jam_f(3))*t/fs);
Jam_4=A_jam*cos(2*pi*(fm+Jam_f(4))*t/fs);  %四个干扰信号

en=3;
if en==1                                   %en为干扰使能
    Jam=Jam_1;                             %单频干扰
elseif en==2
    Jam=Jam_1+Jam_2;
elseif en==3;
    Jam=Jam_1+Jam_2+Jam_3;
elseif en==4;
    Jam=Jam_1+Jam_2+Jam_3+Jam_4;           %多音干扰，模拟窄带干扰
elseif en==0;
    Jam=zeros(1,L);                        %使能为0是无干扰
end 
xif_nsj=xif_sn+Jam;                        %接收到的带干扰和噪声的信号

%频域干扰抑制算法要求对输入信号做fft变换，为了减少频谱泄露，减少干扰信号旁瓣对导航信号产生的干扰，需要fft之前用窗函数进行平滑处理
%窗函数：汉宁窗hanning，海明窗hamming，布拉克曼窗blackman，恺撒窗kaiser，通过w取值选取采用的窗函数
w=4;
if w==0
    window=ones(1,L);                     %不加窗
elseif w==1
    window=hanning(1,L);
elseif w==2
    window=hamming(1,L);
elseif w==3
    window=blackman(1,L);
elseif w==4
    window=blackman(1,L);                 %L点布拉克曼窗
elseif w==5
    window=kaiser(1,L);
end

xif_nsj_win=xif_nsj.*window;              %对接收信号窗函数平滑处理

%干扰抑制算法：固定门限——K值法
%K值法是干扰抑制算法中固定门限法的一种，从刚才的幅频图中可以看出淹没在高斯噪声下的导航信号幅值均大约在80db以下，所以选定固定门限80db作为检测门限，大于80db的频点处直接置零处理
xif_K_fft=fft(xif_nsj_win);
f=(0:L-1)*L/fs;
xif_fft_abs=20*log10(abs(xif_K_fft));
figure(4);
subplot(211);
plot(f,xif_fft_abs);                      %绘制幅频图
xlabel('f/MHz');
ylabel('xif_fft_abs');
Th_K=80;
Jam_index=find(xif_fft_abs>Th_K);         %找出接收信号幅值大于80dB的频点
xif_K_fft(Jam_index)=0;                   %置零
xif_fft_K_abs=20*log10(abs(xif_K_fft));
subplot(212);
plot(f,xif_fft_K_abs);                    %绘制幅频图
xlabel('f/MHz');
ylabel('xif_fft_K_abs');
title('K值法消除单音或窄带干扰');
xif_K_ifft=ifft(xif_K_fft);               %K值法抗干扰处理后ifft输出

%干扰抑制算法：固定门限——一阶矩法
%一阶矩法的门限为Th=q*u，其中u为FFT后均值，q为优化系数，采用1024点的FFT并通过加窗处理减少频谱泄露
xif_one_fft=fft(xif_nsj_win);%L点fft
xif_one_fft_abs=abs(xif_one_fft);
xif_one_fft_absmod=20*log10(xif_one_fft_abs);
figure(5);
subplot(211);
plot(f,xif_one_fft_absmod);               %绘制幅频图
xlabel('f/MHz');
ylabel('xif_one_fft_abs');
u=mean(xif_one_fft_abs,2);                %求L点后信号幅值的期望，mean函数为求期望函数
q=1;                                      %优化系数需试探性取值
Th_one=q*u;                               %检测门限
fft_one_index=find(xif_one_fft_abs>Th_one);%找出大于门限的频点即为干扰所在位置
xif_one_fft(fft_one_index)=0;              %置零
xif_one_fft_abs=20*log10(abs(xif_one_fft));
subplot(212);
plot(f,xif_one_fft_abs);                   %绘制幅频图
xlabel('f/MHz');
ylabel('xif_one_fft_abs');
title('一阶矩法消除单音或窄带干扰');
xif_one_ifft=ifft(xif_one_fft);            %一阶矩抗干扰处理后ifft输出

%干扰抑制算法：固定门限法：中值门限法
%中值门限法也是固定门限法的一种，门限取值为接收信号经过fft后的中间值，即Th=r*median{X[k]},注意不是平均值，而是中间值，r为优化系数
xif_median_fft=fft(xif_nsj_win);           %L点fft
xif_median_fft_abs=abs(xif_median_fft);
xif_median_fft_absmod=20*log10(xif_median_fft_abs);
figure(6);
subplot(211);
plot(f,xif_median_fft_absmod);             %绘制幅频图
xlabel('f/MHz');
ylabel('xif_median_fft_abs');
xif_median_fftsort=sort(xif_median_fft_abs);
median=0.5*(xif_median_fftsort(L/2)+xif_median_fftsort(L/2+1));%将fft变换后的数据从小到达排序并找出中间值
r=4;                                       %优化系数试探性设为1.5
Th_median=r*median;                        %检测门限
fft_median_index=find(xif_median_fft_abs>Th_median);%找出大于门限的频点即为干扰所在位置
xif_median_fft(fft_median_index)=0;        %置零
xif_median_fft_abs=20*log10(abs(xif_median_fft));
subplot(212);
plot(f,xif_median_fft_abs);                %绘制幅频图
xlabel('f/MHz');
ylabel('xif_median_fft_abs');
title('中值门限法消除单音或窄带干扰');
xif_median_ifft=ifft(xif_median_fft);      %中值门限法抗干扰处理后ifft输出

%频域自适门限法
%接收到的导航信号近似为窄带高斯信号，其包络服从瑞丽分布，包络分平方服从指数分布
xif_adjust_fft=fft(xif_nsj_win);%L点fft
xif_adjust_fft_abs=abs(xif_adjust_fft);
xif_adjust_fft_absmod=20*log10(xif_adjust_fft_abs);
figure(7);
subplot(211);
plot(f,xif_adjust_fft_absmod);              %绘制幅频图
xlabel('f/MHz');
ylabel('xif_adjust_fft_abs');
xif_adjust_fft_abs2=abs(xif_adjust_fft).*abs(xif_adjust_fft);
u_adjust=mean(xif_adjust_fft_abs2);        %求估计的指数分布的均值
Th_adjust=4*u_adjust;                      %自适应门限
%根据包络平方的取值范围确定fft的取值
fft_adjust_index1=find(0*Th_adjust<xif_adjust_fft_abs2<Th_adjust);
fft_adjust_index2=find(Th_adjust<xif_adjust_fft_abs2<8*Th_adjust);
fft_adjust_index3=find(8*Th_adjust<xif_adjust_fft_abs2<64*Th_adjust);
fft_adjust_index4=find(64*Th_adjust<xif_adjust_fft_abs2<512*Th_adjust);
fft_adjust_index5=find(xif_adjust_fft_abs2>512*Th_adjust);
xif_adjust_fft(fft_adjust_index1)=xif_adjust_fft(fft_adjust_index1);
xif_adjust_fft(fft_adjust_index2)=xif_adjust_fft(fft_adjust_index2)/8;
xif_adjust_fft(fft_adjust_index3)=xif_adjust_fft(fft_adjust_index3)/64;
xif_adjust_fft(fft_adjust_index4)=xif_adjust_fft(fft_adjust_index4)/512;
xif_adjust_fft(fft_adjust_index5)=0;       %其余值置为0
xif_adjust_fft_abs=20*log10(abs(xif_adjust_fft));
subplot(212);
plot(f,xif_adjust_fft_abs);                %绘制幅频图
xlabel('f/MHz');
ylabel('xif_adjust_fft_abs');
title('自适应门限法消除单音或窄带干扰');
xif_adjust_ifft=ifft(xif_adjust_fft);      %自适应门限法抗干扰处理后ifft输出

