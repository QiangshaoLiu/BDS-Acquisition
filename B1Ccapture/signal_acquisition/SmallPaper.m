%%本程序用于小论文发表，实现了基于ASPeCT的并行码相位捕获算法
%时间：2018年3月
%作者：LSQ

clc;
close all;
%%产生MBOC模拟中频信号，仅包括数据通道
%产生weil码

legendrelength = 10243;
WeilCodelength = 10230;
legendre = zeros(1, legendrelength);
WeilCode = zeros(1, WeilCodelength);
WeilCode_p = zeros(1, WeilCodelength);
phase_diff = 2678;
intercept_pot = 699;
phase_diff_p = 796;   %PRN=1的导频分量主码相位差
intercept_pot_p = 7575;  %PRN=1的导频分量主码截取点

for k = 1 : legendrelength-1 
    for x = 1 : (legendrelength-1)/2
        if mod(k,legendrelength) == mod(x^2, legendrelength)
            legendre(k) = 1;
        end
    end
end

legendre = [0 legendre(1:legendrelength-1)];   %生成legendre序列 

%生成B1C数据分量主码
for k = 0 : WeilCodelength-1     
    WeilCode(k+1) = mod(sum([legendre(mod((k+intercept_pot-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff+intercept_pot-1), legendrelength)+1)]),2);
end

%将Weil码变换成极性码，0表示高电平‘+1’，1表示低电平‘-1’
for k = 1 : WeilCodelength
    if WeilCode(k) == 0
       WeilCode(k) = 1;
    else
       WeilCode(k) = -1;
    end
end

sub_legendrelength = 3607;    %导频分量子码
sub_WeilCodelength = 1800;
sub_legendre = zeros(1, sub_legendrelength);
sub_WeilCode = zeros(1, sub_WeilCodelength);
sub_phase_diff_p = 269;
sub_intercept_pot_p = 1889;

for k = 1 : sub_legendrelength-1 
    for x = 1 : (sub_legendrelength-1)/2
        if mod(k,sub_legendrelength) == mod(x^2, sub_legendrelength)
            sub_legendre(k) = 1;
        end
    end
end

sub_legendre = [0 sub_legendre(1:sub_legendrelength-1)]; %生成子码legendre序列

%生成B1C导频分量子码
for k = 0 : sub_WeilCodelength-1     
    sub_WeilCode(k+1) = mod(sum([sub_legendre(mod((k+sub_intercept_pot_p-1), sub_legendrelength)+1),... 
        sub_legendre(mod((k+sub_phase_diff_p+sub_intercept_pot_p-1), sub_legendrelength)+1)]),2);
end

%生成B1C导频分量的复合码
%不需要，因为仿真时间为10ms，一个主码周期

%生成B1C导频分量主码
for k = 0 : WeilCodelength-1     
    WeilCode_p(k+1) = mod(sum([legendre(mod((k+intercept_pot_p-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff_p+intercept_pot_p-1), legendrelength)+1)]),2);
end

%将Weil码变换成极性码，0表示高电平‘+1’，1表示低电平‘-1’
for k = 1 : WeilCodelength
    if WeilCode_p(k) == 0
       WeilCode_p(k) = 1;
    else
       WeilCode_p(k) = -1;
    end
end

prn_fft = fft(WeilCode,WeilCodelength);
prn_pow = prn_fft.*conj(prn_fft);
prn_cor = ifft(prn_pow);
figure(1)
plot(fftshift(prn_cor));title('Weil码自相关函数');

f_sample = 62*1.023e6;             %采样频率
f_sc_a = 1.023e6 ;           %数据分量子载波速率
f_sc_b = 6*1.023e6 ;         %导频分量子载波速率
T_process = 4e-3;                 %处理时间
T_int = 2e-3;                  %积分时间
Non_Coh_Sums = 5;              %5ms非相干积分
t = 0 : 1/f_sample : T_process - 1/f_sample;

j=sqrt(-1);

Subcarr1 = sign(sin(2*pi*f_sc_a*t));
Subcarr2 = sign(sin(2*pi*f_sc_b*t));
figure(2)
subplot(2,1,1);plot(Subcarr1);title('副载波BOC(1,1)时域波形');
axis([0 f_sample*T_process -1.2 1.2]);
subplot(2,1,2);plot(Subcarr2);title('副载波BOC(6,1)时域波形');
axis([0 f_sample*T_process -1.2 1.2]);

Rc = 1.023e6;       %主码码速率
CodeSample = WeilCode(mod(floor(t*Rc),10230)+1);
CodeSample_p = WeilCode_p(mod(floor(t*Rc),10230)+1);
Subcarr_data = 1/2*CodeSample.*Subcarr1;      
Subcarr_poilt = sqrt(1/11)*CodeSample_p.*Subcarr2 + ...
    j*sqrt(29/44)*CodeSample_p.*Subcarr1;

Subcarr_data_fft = fft(Subcarr_data);
Subcarr_data_pow = Subcarr_data_fft.*conj(Subcarr_data_fft);
Subcarr_data_cor = ifft(Subcarr_data_pow);
Subcarr_data_cor_max = max(Subcarr_data_cor);

figure(3)
N = length(t);
k = 0:N-1;
plot(k, fftshift(Subcarr_data_cor)/Subcarr_data_cor_max);
title('经BOC调制后的导航信号归一化自相关函数');
num_sample = floor(f_sample/Rc);
x_l = N/2-num_sample;
x_h = N/2+num_sample;
axis([x_l x_h -0.6 1.2]);

num = length(Subcarr_data);
delay = 232*num_sample;       %给数据通道信号设定延时
Subcarr_data_delay = [Subcarr_data(delay : num) Subcarr_data(1 : delay-1)];
Subcarr_data_snr = awgn(Subcarr_data_delay, -20);

delay_p = 356*num_sample;
Subcarr_poilt_delay = [Subcarr_poilt(delay_p : num) Subcarr_poilt(1 : delay_p-1)];
Subcarr_poilt_snr = awgn(Subcarr_poilt_delay, -20);

IF = 24.58e6;     %中频频率
fd = 1500;        %多普勒频移
data = Subcarr_data_snr.*sin(2*pi*(IF+fd)*t);    %模拟中频信号
data_row = Subcarr_data_delay.*sin(2*pi*(IF+fd)*t);
data_pow = abs(fftshift(fft(data)));
data_row_pow = abs(fftshift(fft(data_row)));

poilt = Subcarr_poilt_snr.*sin(2*pi*(IF+fd)*t);
poilt_row = Subcarr_poilt_delay.*sin(2*pi*(IF+fd)*t);
poilt_pow = abs(fftshift(fft(poilt)));
poilt_row_pow = abs(fftshift(fft(poilt_row)));

signal_total = Subcarr_data + Subcarr_poilt;
delay_all = 325*num_sample; 
signal_total_delay = [signal_total(delay : num) signal_total(1 : delay-1)];
signal_total_m = signal_total_delay.*sin(2*pi*(IF+fd)*t); 
signal_total_snr = awgn(signal_total_m, -20);

signal_pow = abs(fftshift(fft(signal_total_snr)));
signal_row_pow = abs(fftshift(fft(signal_total_m)));

figure(4)
subplot(3,2,1),plot(data_row_pow);
xlabel('f/Hz');ylabel('Amplitude');
axis([0,3e5,0,3000]),title('Spectrum of data channel signal');
subplot(3,2,2),plot(data_pow);
xlabel('f/Hz');ylabel('Amplitude');
title('With WGN');
subplot(3,2,3),plot(poilt_row_pow);
xlabel('f/Hz');ylabel('Amplitude');
axis([0,3e5,0,3000]),title('Spectrum of pilot channel signal');
subplot(3,2,4),plot(poilt_pow);
xlabel('f/Hz');ylabel('Amplitude');
title('With WGN');

subplot(3,2,5),plot(signal_row_pow);
xlabel('f/Hz');ylabel('Amplitude');
axis([0,3e5,0,3000]),title('Spectrum of navigation signal');
subplot(3,2,6),plot(signal_pow);
xlabel('f/Hz');ylabel('Amplitude');
title('With WGN');

acqSearchStep = 250;      %[Hz]
DopplerRange = 5000;      %[Hz]

FD_vect= -DopplerRange:acqSearchStep:DopplerRange;

% FFT acquisition
idx = 1;
n=0:N-1;
ind_cod = mod(floor(n*Rc/f_sample),WeilCodelength)+1;
SigLOC_tot = WeilCode(ind_cod);
SigLOC_tot_p = WeilCode_p(ind_cod);

C = zeros(length(FD_vect),N);
C_p = zeros(length(FD_vect),N);

for ind_FD= 1:length(FD_vect)
    FD = FD_vect(ind_FD);

    corr = zeros(1,N) + j*zeros(1,N); 

    m= 0:N-1;
    argx = 2*pi*(IF+FD)/f_sample;
    carrI = cos(argx*m);
    carrQ = sin(argx*m);
      
    SigLOC = SigLOC_tot(1:N);
    SigLOCFFT = conj(fft(SigLOC,N)); 
  
    BocLOCFFT = conj(fft(Subcarr1,N));

    SigIN = data(1:N);
    SigOUTI = SigIN .* carrI;
    SigOUTQ = SigIN .* carrQ;
       
    corr = corr - abs(ifft(fft(SigOUTI,N).*(BocLOCFFT))).^2 + ...
            abs(ifft(fft(SigOUTI,N).*(SigLOCFFT))).^2 - ...
            abs(ifft(fft(SigOUTQ,N).*(BocLOCFFT))).^2 + ...
            abs(ifft(fft(SigOUTQ,N).*(SigLOCFFT))).^2;   
            
    C(idx,:) = corr;
    idx = idx+1;
end

index = 1;
for ind_FD= 1:length(FD_vect)
    FD = FD_vect(ind_FD);

    corr = zeros(1,N) + j*zeros(1,N); 

    m= 0:N-1;
    argx = 2*pi*(IF+FD)/f_sample;
    carrI = cos(argx*m);
    carrQ = sin(argx*m);
      
    SigLOC_p = SigLOC_tot_p(1:N);
    SigLOCFFT_p = conj(fft(SigLOC_p,N)); 
  
    BocLOCFFT_p = conj(fft(Subcarr1,N));
    BocLOCFFT_p2 = conj(fft(j*Subcarr2,N));

    SigIN_p = poilt(1:N);
    SigOUTI_p = SigIN_p .* carrI;
    SigOUTQ_p = SigIN_p .* carrQ;
       
    corr = corr - abs(ifft(fft(SigOUTI_p,N).*(BocLOCFFT_p))).^2 + ...
            abs(ifft(fft(SigOUTI_p,N).*(SigLOCFFT_p))).^2 - ...
            abs(ifft(fft(SigOUTQ_p,N).*(BocLOCFFT_p))).^2 + ...
            abs(ifft(fft(SigOUTQ_p,N).*(SigLOCFFT_p))).^2 - ...
        abs(ifft(fft(SigOUTI_p,N).*(BocLOCFFT_p2))).^2 - ...
        abs(ifft(fft(SigOUTQ_p,N).*(BocLOCFFT_p2))).^2;   %QMBOC
            
    C_p(index,:) = corr;
    index = index+1;
end

%Find the main peak in the correlation floor and the corresponding frequency bin index
[bb, ind_mixf] = max(max(C'));
[bb, ind_mixc] = max(max(C));
[pp, ind_mixf_p] = max(max(C_p'));
[pp, ind_mixc_p] = max(max(C_p));

%code_phase = floor((int_n - ind_mixc)/num_sample);
code_phase = ceil((N - ind_mixc)/num_sample);
doppler =(ind_mixf-1)*acqSearchStep - DopplerRange;   %[HZ]

code_phase_p = ceil((N - ind_mixc_p)/num_sample);
doppler_p =(ind_mixf_p-1)*acqSearchStep - DopplerRange;   %[HZ]

%画三维图
[C_y, C_x]=size(C);    %C_y为码相位，C_x为多普勒频移
X=1:C_x;Y=1:C_y;       %X为码相位，Y为多普勒频移
[x,y]=meshgrid(X,Y);
figure(5)
mesh((N-x)/num_sample,(y-1)*acqSearchStep - DopplerRange,C);  %坐标变换
hold on
plot3(code_phase,doppler,bb,'k.','markersize',20);     %小黑点标记
text(code_phase,doppler,bb,['X:',num2str(code_phase),char(10), ...
    'Y:',num2str(doppler),char(10),'Z:',num2str(bb),char(10)]);  %标值
xlabel('Code phase/chip');ylabel('Doppler shift/Hz');zlabel('Correlation value');
title('The acqusition result of the Beidou B1C data channel signal');


[C_y_p, C_x_p]=size(C_p);    %C_y为码相位，C_x为多普勒频移
X_p=1:C_x_p;Y_p=1:C_y_p;       %X为码相位，Y为多普勒频移
[x_p,y_p]=meshgrid(X_p,Y_p);
figure(6)
mesh((N-x_p)/num_sample,(y_p-1)*acqSearchStep - DopplerRange,C_p);  %坐标变换
hold on
plot3(code_phase_p,doppler_p,pp,'k.','markersize',20);     %小黑点标记
text(code_phase_p,doppler_p,pp,['X:',num2str(code_phase_p),char(10), ...
    'Y:',num2str(doppler_p),char(10),'Z:',num2str(pp),char(10)]);  %标值
xlabel('Code phase/chip');ylabel('Doppler shift/Hz');zlabel('Correlation value');
title('The acqusition result of the Beidou B1C pilot channel signal');


