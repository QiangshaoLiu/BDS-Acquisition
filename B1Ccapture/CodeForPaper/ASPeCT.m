%Author:LSQ
%Date:2019/4
%Description: ASPeCT捕获算法验证.

clc;
close all;

set(0,'defaultfigurecolor','w'); %将仿真图背景设置为白色

%对BOC(1,1)和QMBOC(6,1,4/33)的边锋消除效果进行了分析
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

f_sample = 36*1.023e6;             %采样频率
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

n=0:length(t)-1;
ind_cod = mod(floor(n*Rc/f_sample),WeilCodelength)+1;
SigLOC_tot_d = WeilCode_d(ind_cod);
SigLOC_tot_p = WeilCode_p(ind_cod);

%AsPect方法所得到的自相关函数
[BOCPRN_corr lag] = xcorr(Subcarr1.*SigLOC_tot_d, SigLOC_tot_d, 'coeff');  %BOC信号与PRN码的互相关函数
BOC_corr = xcorr(Subcarr1.*SigLOC_tot_d, 'coeff');    %BOC的自相关函数
Aspect = abs(BOC_corr).^2 - 0.8*abs(BOCPRN_corr).^2;
Aspect_1 = abs(BOC_corr).^2 - 0.9*abs(BOCPRN_corr).^2;
Aspect_2 = abs(BOC_corr).^2 - 1*abs(BOCPRN_corr).^2;
Aspect_3 = abs(BOC_corr).^2 - 1.1*abs(BOCPRN_corr).^2;
Aspect_4 = abs(BOC_corr).^2 - 1.2*abs(BOCPRN_corr).^2;
Aspect_5 = abs(BOC_corr).^2 - 1.3*abs(BOCPRN_corr).^2;

%横坐标尺度变换
index = lag / num_sample;

figure(1)
plot(index, Aspect, 'r', index, BOCPRN_corr, 'g', index, BOC_corr, 'b');
xlabel('相位延时(码片)');
ylabel('归一化相关函数');
axis([-2.5 2.5 -0.6 1.4]);
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
legend('重构的相关函数','BOC-PRN互相关函数','BOC(1,1)自相关函数');

%%用于证明利用ASPeCT捕获算法对QMBOC信号进行处理会存在自相关边峰消除不干净的情况
QMBOCPRN_corr = xcorr(B1C_poilt, SigLOC_tot_p, 'coeff'); 
QMBOC_corr = xcorr(B1C_poilt, 'coeff');
Aspect_qmboc = abs(QMBOC_corr).^2 - abs(QMBOCPRN_corr).^2;
figure(2)
plot(index, Aspect_qmboc, 'r', index, QMBOCPRN_corr, 'g', index, QMBOC_corr, 'b');
xlabel('相位延时(码片)');
ylabel('归一化相关函数');
axis([-2.5 2.5 -0.6 1.4]);
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
legend('重构的相关函数','QMBOC-PRN互相关函数','QMBOC(6,1,4/33)自相关函数');
title('QMBOC(6,1,4/33)的重构相关函数');
