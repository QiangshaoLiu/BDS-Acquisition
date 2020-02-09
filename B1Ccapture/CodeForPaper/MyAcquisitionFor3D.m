%Author:LSQ
%Date:2019/4
%Description: 为学位论文设计的B1C信号捕获算法，基于PCF方法进行了改进，给出了捕获结果3D图.

clc;
close all;

set(0,'defaultfigurecolor','w'); %将仿真图背景设置为白色

%%本程序对北斗B1C信号的导频和数据分量进行联合捕获
%仿真参数设置
f_sample = 36*1.023e6;             %采样频率
f_sc_a = 1.023e6 ;                 %数据分量子载波速率
f_sc_b = 6*1.023e6 ;               %导频分量子载波速率
Rc = 1.023e6;                      %主码码速率
T_process = 20e-3;                 %处理时间
T_int = 10e-3;                      %相关运算时间
Non_Coh_Sums = 2;                  %(Non_Coh_Sums*T_int)ms非相干积分时间
t = 0 : 1/f_sample : T_process - 1/f_sample;
n = 0:length(t)-1;                 
j=sqrt(-1);
pi = 3.141592654;                  %圆周率
Num_int = floor(f_sample * T_int); %相干积分时间所对应的采样点数

%%模拟产生接收信号
subcarr1 = sign(sin(2*pi*f_sc_a*t));
subcarr1(1) = 1;
subcarr2 = sign(sin(2*pi*f_sc_b*t));
subcarr2(1) = 1;
code_r = generatecode(2);           %接收信号由PRN=2的扩频码序列调制
code_data_r = generatedatacode(2);  %接收信号的数据分量
codeSample_r = code_r(mod(floor(t*Rc),10230)+1);
codeSample_data_r = code_data_r(mod(floor(t*Rc),10230)+1);
Qmboc_p = sqrt(1/11)*codeSample_r.*subcarr2 + ...
    j*sqrt(29/44)*codeSample_r.*subcarr1;
Boc_d = 1/2*codeSample_data_r.*subcarr1;    %默认处理时间内的导航电文全为1

Signal_B1C = Qmboc_p + Boc_d;

code_sample = floor(f_sample/Rc);   %单个码片所对应的采样数
num_signal = length(Signal_B1C);
delay = 306*code_sample;            %给伪码设定码相位延时
Signal_B1C_delay = [Signal_B1C(delay : num_signal) Signal_B1C(1 : delay-1)];

IF = 24.58e6;     %中频频率
fd = 1200;        %多普勒频移
signal_IF = Signal_B1C_delay.*cos(2*pi*(IF+fd)*t); %模拟中频信号，只考虑IQ分量的I分量

signal = awgn(signal_IF, -37);    %加高斯白噪声

%%基于PCF的北斗B1C信号捕获算法
FdSearchStep = 50;      %[Hz]
DopplerRange = 5000;      %[Hz]

FdVect= -DopplerRange:FdSearchStep:DopplerRange;     %多普勒频移搜索范围

%产生本地测距码序列
prn_p = generatecode(2);
prn_d = generatedatacode(2);
index_code = mod(floor(Rc*t),10230)+1;
prn_local = prn_p(index_code);
prn_d_local = prn_d(index_code);

%%导频信号QMBOC(6,1,4/33)的BOC(1,1)
idx1 = mod(floor(12*Rc*t),12)+1;
prn1_qmboc11 = [j*sqrt(2),j*sqrt(2),j*sqrt(2),j*sqrt(2),j*sqrt(2),...
    j*sqrt(2),0,0,0,0,0,0];
s1_qmboc11 = prn1_qmboc11(idx1).*prn_local;
[g1_qmboc11 x]= xcorr(Qmboc_p, s1_qmboc11, 'coeff');
prn12_qmboc11 = [0,0,0,0,0,0,j*sqrt(2),j*sqrt(2),j*sqrt(2),j*sqrt(2),...
    j*sqrt(2),j*sqrt(2)];
s12_qmboc11 = prn12_qmboc11(idx1).*prn_local;
%g12_qmboc11 = xcorr(Qmboc_p, s12_qmboc11, 'coeff');

%corr_sum_qmboc11 = abs(g1_qmboc11)+abs(g12_qmboc11)-abs(g1_qmboc11+g12_qmboc11);
%corr_qmboc11 = xcorr(Qmboc_p, Qmboc_p, 'coeff');

%%数据分量BOC(1,1)以及导频分量BOC(6,1)
prn1_qmboc61 = [sqrt(6),0,0,0,0,0,0,0,0,0,0,0];
s1_qmboc61 = prn1_qmboc61(idx1).*prn_local;
%[g1_qmboc61 x]= xcorr(Qmboc_p, s1_qmboc61, 'coeff');
prn12_qmboc61 = [0,0,0,0,0,0,0,0,0,0,0,sqrt(6)];
s12_qmboc61 = prn12_qmboc61(idx1).*prn_local;
%g12_qmboc61 = xcorr(Qmboc_p, s12_qmboc61, 'coeff');

%corr_sum_qmboc61 = abs(g1_qmboc61)+abs(g12_qmboc61)-abs(g1_qmboc61+g12_qmboc61);
%corr_qmboc61 = xcorr(Qmboc_p, Qmboc_p, 'coeff');

%corr_sum_qmboc = corr_sum_qmboc11 + corr_sum_qmboc61;

prn1_boc11 = [1,1,1,1,1,1,0,0,0,0,0,0];
s1_boc11 = prn1_boc11(idx1).*prn_d_local;
%[g1_boc11 x]= xcorr(Boc_d, s1_boc11, 'coeff');
prn2_boc11 = [0,0,0,0,0,0,1,1,1,1,1,1];
s2_boc11 = prn2_boc11(idx1).*prn_d_local;
%g2_boc11 = xcorr(Boc_d, s2_boc11, 'coeff');

%corr_sum_boc11 = abs(g1_boc11)+abs(g2_boc11)-abs(g1_boc11+g2_boc11);
%corr_boc11 = xcorr(Boc_d, Boc_d, 'coeff');

%横坐标尺度变换
index = x / floor(f_sample/Rc);

% %作图导频分量
% figure(1)
% plot(index, corr_qmboc11,'b',index,corr_sum_qmboc11,'m');
% hold on;
% plot(index, corr_sum_qmboc,'r');
% xlabel('码片');
% ylabel('归一化自相关函数');
% legend('ACF','PCF','Pilot');
% axis([-1.5 1.5 -0.5 2]);
% 
% %作图数据分量
% figure(2)
% plot(index, corr_boc11,'b',index,corr_sum_boc11,'m');
% xlabel('码片');
% ylabel('归一化自相关函数');
% legend('ACF','PCF');
% axis([-1.5 1.5 -0.5 1.5]);

%%以下为PCF捕获算法验证
%生成矩阵用于存相关结果
C = zeros(length(FdVect),Num_int);     %用于所有码片的相关结果
C_t = zeros(1,Num_int);     %用于存门限值
idx = 1;     %矩阵行数

for ind_FD= 1:length(FdVect)
    corr_temp = zeros(1,Num_int) ;
    fd_ind = FdVect(ind_FD);
    %本地载波
    m = 1:Num_int;
    carrI = cos(2*pi*(IF+fd_ind)*m/f_sample);
    carrQ = sin(2*pi*(IF+fd_ind)*m/f_sample);
    for M = 0 : (Non_Coh_Sums - 1)
    %下变频
    SigIN = signal(1+M*Num_int : Num_int+M*Num_int);
    SigOUTI = SigIN .* carrI;
    SigOUTQ = SigIN .* carrQ;
    %本地导频码
    S1_qmboc11 = s1_qmboc11(1:Num_int);
    S12_qmboc11 = s12_qmboc11(1:Num_int);
    S1_qmboc61 = s1_qmboc61(1:Num_int);
    S12_qmboc61 = s12_qmboc61(1:Num_int);
    
    PRNLOCFFT_boc11_E = conj(fft(S1_qmboc11));
    PRNLOCFFT_boc11_L = conj(fft(S12_qmboc11));
    PRNLOCFFT_boc61_E = conj(fft(S1_qmboc61));
    PRNLOCFFT_boc61_L = conj(fft(S12_qmboc61));
    %本地数据码
    S1_boc11 = s1_boc11(1:Num_int);
    S2_boc11 = s2_boc11(1:Num_int);
    
    PRNLOCFFT_boc11_d_E = conj(fft(S1_boc11));
    PRNLOCFFT_boc11_d_L = conj(fft(S2_boc11));
    
    %对基带信号进行FFT
    SigOUT = SigOUTI + SigOUTQ;
    Signal_fft = fft(SigOUT);
    
    %重构相关函数
    R_boc_prn_E_11 = Signal_fft.*PRNLOCFFT_boc11_E;
    R_boc_prn_L_11 = Signal_fft.*PRNLOCFFT_boc11_L;
    R_boc_prn_E_61 = Signal_fft.*PRNLOCFFT_boc61_E;
    R_boc_prn_L_61 = Signal_fft.*PRNLOCFFT_boc61_L;
    R_boc_prn_d_E_11 = Signal_fft.*PRNLOCFFT_boc11_d_E;
    R_boc_prn_d_L_11 = Signal_fft.*PRNLOCFFT_boc11_d_L;
    
    R_E_11 = ifft(R_boc_prn_E_11);
    R_L_11 = ifft(R_boc_prn_L_11);
    R_E_61 = ifft(R_boc_prn_E_61);
    R_L_61 = ifft(R_boc_prn_L_61);
    R_d_E_11 = ifft(R_boc_prn_d_E_11);
    R_d_L_11 = ifft(R_boc_prn_d_L_11);

    R_EL_11 = R_E_11 + R_L_11;
    R_EL_61 = R_E_61 + R_L_61;
    R_d_EL_11 = R_d_E_11 + R_d_L_11;
    
    corr_temp = corr_temp + abs(R_E_11) + abs(R_L_11) - abs(R_EL_11)...
        + abs(R_E_61) + abs(R_L_61) - abs(R_EL_61)...
        + abs(R_d_E_11) + abs(R_d_L_11) - abs(R_d_EL_11);
    end
       
    C(idx,:) = corr_temp;
    idx = idx + 1;
end

[value, ind_mixf] = max(max(C'));
[value, ind_mixc] = max(max(C));

code_phase = (Num_int-ind_mixc)/code_sample;
doppler =(ind_mixf-1)*FdSearchStep - DopplerRange;   %[HZ]

%自适应捕获检测判决
Num_code = 12;                %被检测码片单元周围的参考码片单元数目
ThresholdFactor = 9.34;      %虚警率为10e-6的门限比例因子
Z = 0;                        %功率估计值
for i=1:6
    Z = Z + C(ind_mixf,ind_mixc+i*code_sample)+C(ind_mixf,ind_mixc-i*code_sample);    
end
Z_aver = Z/Num_code;

Threshold = ThresholdFactor*Z_aver - 3e5;    %得到自适应门限值

for i=1:Num_int
     C_t(i) = Threshold;
end

if C(ind_mixf,ind_mixc) > Threshold
    data = sprintf('The acquisition result\n Code phase:%f 码片\nDoppler frequency:%f Hz\nThreshold:%f \n',...
        code_phase,doppler,Threshold);
    disp(data);
else
    data = sprintf('Acquisition failed!\n');
    disp(data);
end

%画三维图
[C_y, C_x]=size(C);    %C_x为码相位，C_y为多普勒频移
X=1:C_x;Y=1:C_y;       %X为码相位，Y为多普勒频移
[x,y]=meshgrid(X,Y);
% figure(3)
% mesh((Num_int-x)/code_sample,(y-1)*FdSearchStep - DopplerRange,C);  %坐标变换
% hold on;
% plot3(code_phase,doppler,value,'k.','markersize',20);     %小黑点标记
% text(code_phase,doppler,value,['X:',num2str(code_phase,3),char(10), ...
%      'Y:',num2str(doppler,4),char(10),'Z:',num2str(value),char(10)]);  %标值
% xlabel('码相位延时(码片)');ylabel('多普勒频移(Hz)');zlabel('相关值');
% axis([0 10230 -5000 5000 0 value+1e4]);
C_x_2D = (doppler+DopplerRange)/FdSearchStep+1;
C_2D = C(C_x_2D,:);
figure(3)
plot((Num_int-x)/code_sample, C_2D,'b',(Num_int-x)/code_sample, C_t,'r--');  %画出门限值
legend('多普勒频移：1200Hz');
%axis([code_phase-2 code_phase+2 0 value+1e4]);
axis([0 10230 0 value+1e4]);
xlabel('码相位延时(码片)');ylabel('相关值');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

