%本程序为毕业设计捕获算法最终方案的仿真代码，利用MATLAB模拟产生数字中频信号，
%最终得到捕获结果的三维图，二维图以及相关信号的自相干函数
%作者：LSQ
%时间：2018年4月

clc;
close all;

%产生Weil码，PRN号为1至10,数据分量和导频分量主码
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
        intercept_pot_p = 7575;    %超出范围，默认PRN=1
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

legendre = [0 legendre(1:legendrelength-1)];   %生成legendre序列 

%生成B1C数据/导频分量主码
for k = 0 : WeilCodelength-1     
    WeilCode_d(k+1) = mod(sum([legendre(mod((k+intercept_pot_d-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff_d+intercept_pot_d-1), legendrelength)+1)]),2);
    WeilCode_p(k+1) = mod(sum([legendre(mod((k+intercept_pot_p-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff_p+intercept_pot_p-1), legendrelength)+1)]),2);
end

%将Weil码变换成极性码，0表示高电平‘+1’，1表示低电平‘-1’
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

f_sample = 20*1.023e6;             %采样频率
f_sc_a = 1.023e6 ;                 %数据分量子载波速率
f_sc_b = 6*1.023e6 ;               %导频分量子载波速率
T_process = 15e-3;                 %处理时间
T_int = 2e-3;                      %相关运算时间
Non_Coh_Sums = 5;                  %(Non_Coh_Sums*T_int)ms非相干积分时间
t = 0 : 1/f_sample : T_process - 1/f_sample;

j=sqrt(-1);

%N = length(t);     %不做非相干积分，直接对处理时间内的数据进行ASPeCT处理
N = floor(f_sample * T_int);

Subcarr1 = sign(sin(2*pi*f_sc_a*t));
Subcarr2 = sign(sin(2*pi*f_sc_b*t));

Rc = 1.023e6;       %主码码速率
CodeSample_d = WeilCode_d(mod(floor(t*Rc),10230)+1);
CodeSample_p = WeilCode_p(mod(floor(t*Rc),10230)+1);
B1C_data = 1/2*CodeSample_d.*Subcarr1;      
B1C_poilt = sqrt(1/11)*CodeSample_p.*Subcarr2 + ...
    j*sqrt(29/44)*CodeSample_p.*Subcarr1;

%数据通道和导频通道的组合信号
Signal_B1C = B1C_data + B1C_poilt;

num_sample = floor(f_sample/Rc);
num = length(Signal_B1C);
delay = 273*num_sample;       %给信号设定延时
Signal_B1C_delay = [Signal_B1C(delay : num) Signal_B1C(1 : delay-1)];

IF = 24.58e6;     %中频频率
fd = 1500;        %多普勒频移
Signal_B1C_out = Signal_B1C_delay.*cos(2*pi*(IF+fd)*t); %模拟中频信号，只考虑IQ分量的I分量

data = awgn(Signal_B1C_out, -25);    %加高斯白噪声

acqSearchStep = 250;      %[Hz]
DopplerRange = 5000;      %[Hz]

FD_vect= -DopplerRange:acqSearchStep:DopplerRange;

% FFT acquisition
idx = 1;
n=0:length(t)-1;
ind_cod = mod(floor(n*Rc/f_sample),WeilCodelength)+1;
SigLOC_tot_d = WeilCode_d(ind_cod);
SigLOC_tot_p = WeilCode_p(ind_cod);

%AsPect方法所得到的自相关函数
[BOCPRN_corr lag] = xcorr(Subcarr1.*SigLOC_tot_d, SigLOC_tot_d, 'coeff');  %BOC信号与PRN码的互相关函数
PRN_corr = xcorr(SigLOC_tot_d, 'coeff');    %PRN码的自相关函数
Aspect = abs(PRN_corr).^2 - abs(BOCPRN_corr).^2;

%横坐标尺度变换
index = lag / num_sample;

figure(1)
plot(index, Aspect, 'r', index, BOCPRN_corr, 'g', index, PRN_corr, 'b');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-2.5 2.5 -0.6 1.4]);
legend('基于ASPeCT的相关函数','BOC-PRN互相关函数','MBOC信号自相关函数');
title('基于ASPeCT的归一化相关函数');
   
C = zeros(length(FD_vect),N);
C_non = zeros(length(FD_vect),N);  %不做非相干积分
C_p = zeros(length(FD_vect),N);

for ind_FD= 1:length(FD_vect)
    FD = FD_vect(ind_FD);
    corr = zeros(1,N) + j*zeros(1,N);
    corr_non = zeros(1,N) + j*zeros(1,N);
    corr_p = zeros(1,N) + j*zeros(1,N);
    %本地载波
    m = 0:N-1;
    argx = 2*pi*(IF+FD)/f_sample;
    carrI = cos(argx*m);
    carrQ = sin(argx*m);

    %非相干累加
    for M = 0 : (Non_Coh_Sums - 1);
    %for M = 0 : (1 - 1);
    %下变频
    SigIN = data(N*M + 1 : N*M + N);
    SigOUTI = SigIN .* carrI;
    SigOUTQ = SigIN .* carrQ;
    %本地码和本地副载波
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
    
    %ASPeCT: BOC-BOC 减 BOC-PRN
    corr = corr + abs(ifft(fft(SigOUTI,N).*BOCLOCFFT_conj)).^2 - ...
        abs(ifft(fft(SigOUTI,N).*PRNLOCFFT_conj)).^2 + ...
        abs(ifft(fft(SigOUTQ,N).*BOCLOCFFT_conj)).^2 - ...
        abs(ifft(fft(SigOUTQ,N).*PRNLOCFFT_conj)).^2;
    %不做相干积分
    corr_non = abs(ifft(fft(SigOUTI,N).*BOCLOCFFT_conj)).^2 - ...
        abs(ifft(fft(SigOUTI,N).*PRNLOCFFT_conj)).^2 + ...
        abs(ifft(fft(SigOUTQ,N).*BOCLOCFFT_conj)).^2 - ...
        abs(ifft(fft(SigOUTQ,N).*PRNLOCFFT_conj)).^2; 
    %导频信号捕获
    corr_p = corr_p + abs(ifft(fft(SigOUTI,N).*conj(fft(BOCLOC_p,N)))).^2 - ...
        abs(ifft(fft(SigOUTI,N).*conj(fft(PRNLOC_p,N)))).^2 + ...
        abs(ifft(fft(SigOUTQ,N).*conj(fft(BOCLOC_p,N)))).^2 - ...
        abs(ifft(fft(SigOUTQ,N).*conj(fft(PRNLOC_p,N)))).^2;    
    end
    
    C(idx,:) = corr;
    C_non(idx,:) = corr_non;
    C_p(idx,:) = corr_p;
    idx = idx+1;
end

%Find the main peak in the correlation floor and the corresponding frequency bin index
[bb, ind_mixf] = max(max(C'));
[bb, ind_mixc] = max(max(C));

%code_phase = floor((int_n - ind_mixc)/num_sample);
code_phase = ceil((N - ind_mixc)/num_sample);
doppler =(ind_mixf-1)*acqSearchStep - DopplerRange;   %[HZ]

%画三维图
[C_y, C_x]=size(C);    %C_y为码相位，C_x为多普勒频移
X=1:C_x;Y=1:C_y;       %X为码相位，Y为多普勒频移
[x,y]=meshgrid(X,Y);
figure(2)
mesh((N-x)/num_sample,(y-1)*acqSearchStep - DopplerRange,C);  %坐标变换
hold on
plot3(code_phase,doppler,bb,'k.','markersize',20);     %小黑点标记
text(code_phase,doppler,bb,['X:',num2str(code_phase),char(10), ...
    'Y:',num2str(doppler),char(10),'Z:',num2str(bb),char(10)]);  %标值
xlabel('码相位/码');ylabel('多普勒频移/Hz');zlabel('相关值');
title('北斗B1C信号捕获结果');

%画二维图
C_x_2D = (doppler+DopplerRange)/acqSearchStep+1;
C_2D = C(C_x_2D,:);
C_2D_non = C_non(C_x_2D,:);
C_2D_p = C_p(C_x_2D,:);
figure(3)
fuc1 = plot((N-x)/num_sample, C_2D, 'r');
hold on;
fuc2 = plot((N-x)/num_sample, C_2D_non, 'b');
hold on;
fuc3 = plot((N-x)/num_sample, C_2D_p, 'g');
legend([fuc1(1),fuc2(1),fuc3(1)],'改进的ASPeCT捕获算法','基于ASPeCT捕获算法',...
    '导频通道捕获算法');
axis([code_phase-6 code_phase+6 -3e8 bb]);
xlabel('码相位/码');ylabel('相关值');
title('捕获结果的Y-Z平面');

