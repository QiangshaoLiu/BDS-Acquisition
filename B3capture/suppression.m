%-------------------------------------------------------------------------
% 化腐朽为神奇
% 作者：LSQ
% 日期：2016年11月23日
%----------------------程序说明-------------------------------------------
% 本程序主要用于窄带干扰抑制仿真,使用时域窄带干扰抑制方法
% 使用了有限冲击响应（FIR）陷波和无限冲击响应（IIR）陷波方法
%----------------------基本参数-------------------------------------------
% 采样率：2000Hz 干扰类型：单音干扰 宽带信号：-20db的高斯白噪声
% ---------------------FIR滤波器参数--------------------------------------
% 带阻截止频率：380Hz/620Hz   带通频段：400Hz至600Hz
%-------------------------------------------------------------------------
%时域窄带干扰抑制方法
clear;fs=2000;t=(1:8192)/fs;
x=wgn(1,8192,-20);
x_noise=x+0.01*cos(2*pi*500*t);%加入窄带干扰
L=length(x);N=2^(nextpow2(L));X=fft(x,N);X_NOISE=fft(x_noise,N);
figure(1);subplot(2,1,1);plot((0:N-1)*fs/L,abs(X));
grid on;title('加干扰前信号频谱');xlabel('f');ylabel('振幅|H(e^jw)|');%原始信号
subplot(2,1,2);plot((0:N-1)*fs/L,abs(X_NOISE));%查看信号频谱
grid on;title('加干扰后信号频谱');xlabel('f');ylabel('振幅|H(e^jw)|');
 
%高斯白噪声通过带通滤波器模拟宽带信号
Ap=1;As=60;% 定义通带及阻带衰减
dev=[10^(-As/20),(10^(Ap/20)-1)/(10^(Ap/20)+1),10^(-As/20)];% 计算偏移量
mags=[0,1,0];% 带通
fcuts=[380,400,600,620];% 边界频率
[N,Wn,beta,ftype]=kaiserord(fcuts,mags,dev,fs);% 估算FIR滤波器阶数
hh2=fir1(N,Wn,ftype,kaiser(N+1,beta));% FIR滤波器设计
x_2=filter(hh2,1,x);% 滤波
x_2(1:ceil(N/2))=[];% 群延时N/2，删除无用信号部分
 
x_noise_2=filter(hh2,1,x_noise);% 滤波
x_noise_2(1:ceil(N/2))=[];% 群延时N/2，删除无用信号部分
 
L=length(x_2);N=2^(nextpow2(L));X_2=fft(x_2,N);X_NOISE_2=fft(x_noise_2,N);
figure(2);subplot(3,1,1);plot((0:N-1)*fs/L,abs(X_2));
grid on;title('窄带干扰前的信号');xlabel('f');ylabel('振幅|H(e^jw)|');
axis([250 750 -1 50]);
subplot(3,1,2);plot((0:N-1)*fs/L,abs(X_NOISE_2));% 查看信号频谱
axis([250 750 -1 50]);
grid on;title('窄带干扰后的信号');xlabel('f');ylabel('振幅|H(e^jw)|');
 
%IIR陷波器
Ap2=1;As2=20;
Wp2=[497/1000 503/1000];
Ws2=[499/1000 501/1000];
[N2,Wc2]=buttord(Wp2,Ws2,Ap2,As2);
[b,a]=butter(N2,Ws2,'stop');
bb=filter(b,a,x_noise_2);
bb(1:ceil(N/2))=[];% 群延时N/2，删除无用信号部分
BB=fft(bb,N);
subplot(3,1,3);plot((0:N-1)*fs/L,abs(BB));% 查看信号频谱
axis([250 750 -1 25]);
grid on;title('窄带干扰抑制后的信号');xlabel('f');ylabel('振幅|H(e^jw)|');
