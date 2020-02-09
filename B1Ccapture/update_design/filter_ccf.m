%本程序用于验证filter捕获算法

clc;
close all;

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

f_sample = 40*1.023e6;             %采样频率
f_sc_a = 1.023e6 ;                 %数据分量子载波速率
f_sc_b = 6*1.023e6 ;               %导频分量子载波速率
T_process = 10e-3;                 %处理时间

t = 0 : 1/f_sample : T_process - 1/f_sample;

j=sqrt(-1);

Subcarr1 = sign(sin(2*pi*f_sc_a*t));
Subcarr2 = sign(sin(2*pi*f_sc_b*t));

Rc = 1.023e6;       %主码码速率

% FFT acquisition
idx = 1;
n=0:length(t)-1;
ind_cod = mod(floor(n*Rc/f_sample),WeilCodelength)+1;
SigLOC_tot_d = WeilCode(ind_cod);
SigLOC_tot_p = WeilCode_p(ind_cod);

%AsPect方法所得到的重构相关函数
%数据分量
beta = 1.2;
[BOCPRN_corr lag] = xcorr(Subcarr1.*SigLOC_tot_d, SigLOC_tot_d, 'coeff');  %BOC信号与PRN码的互相关函数
BOC_corr = xcorr(Subcarr1.*SigLOC_tot_d, 'coeff');    %BOC的自相关函数
Aspect = abs(BOC_corr).^2 - abs(BOCPRN_corr).^2;
Aspect_b = abs(BOC_corr).^2 - beta*abs(BOCPRN_corr).^2;
%导频分量
BOCPRN_corr_p = xcorr(-j*sqrt(4/33)*Subcarr2.*SigLOC_tot_p+sqrt(29/33)*Subcarr1.*SigLOC_tot_p, SigLOC_tot_p, 'coeff');  %BOC信号与PRN码的互相关函数
BOC_corr_p = xcorr(-j*sqrt(4/33)*Subcarr2.*SigLOC_tot_p+sqrt(29/33)*Subcarr1.*SigLOC_tot_p, 'coeff');    %BOC的自相关函数
Aspect_p = abs(BOC_corr_p).^2 - abs(BOCPRN_corr_p).^2;
Aspect_p_b = abs(BOC_corr_p).^2 - beta*abs(BOCPRN_corr_p).^2;

num_sample = floor(f_sample/Rc);

%横坐标尺度变换
index = lag / num_sample;

figure(1)
subplot(1,2,1);
plot(index, sqrt(Aspect.^2), 'b', index,sqrt(BOC_corr.^2), 'r');
hold on;
plot(index,sqrt(Aspect_b.^2),'g',index, BOCPRN_corr,'k');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1.5 1.5 -0.5 1.2]);
legend('ASPeCT beta=1.0','BOC(1,1)自相关','ASPeCT beta=1.2','BOC-PRN互相关');
title('数据分量的自相关边峰消除前后对比');
subplot(1,2,2);
plot(index, sqrt(Aspect_p.^2), 'b', index, sqrt(BOC_corr_p.^2), 'r');
hold on;
plot(index,sqrt(Aspect_p_b.^2),'g',index, BOCPRN_corr_p,'k');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1.5 1.5 -0.5 1.2]);
legend('ASPeCT beta=1.0','QMBOC(6,1,4/33)自相关','ASPeCT beta=1.2','QMBOC-PRN互相关');
title('导频分量的自相关边峰消除前后对比');

%用于证明ASPeCT方法仅适用于BOC(n，n)信号
BOCPRN_corr_p_n = xcorr(Subcarr2.*SigLOC_tot_p, SigLOC_tot_p, 'coeff'); %BOC(6,1)与PRN码互相关
BOC_corr_p_n = xcorr(Subcarr2.*SigLOC_tot_p, 'coeff'); 
Aspect_p_n = abs(BOC_corr_p_n).^2 - abs(BOCPRN_corr_p_n).^2;

figure(2)
plot(index, sqrt(Aspect_p_n.^2), 'r', index,sqrt(BOC_corr_p_n.^2), 'b');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1.5 1.5 -0.5 1.2]);
legend('互相关捕获算法','BOC(6,1)自相关');
title('互相关法重构BOC(6,1)信号自相关函数');

%基于filter的方法所得到的重构相关函数
%数据分量
PRN_num =   floor(1000/2000*1/Rc*f_sample);                               %超前/滞后半个码片所对应的采样点数
zero_filling = ones(1, PRN_num-1);
a_num = floor(1150/2000*1/Rc*f_sample); 
a_filling = ones(1, a_num-1);
PRN_lag = [SigLOC_tot_d(PRN_num:end) zero_filling];        %滞后半个码片的PRN
PRN_advance = [a_filling SigLOC_tot_d(1:end-a_num+1)];%超前半个码片的PRN
BOCPRN_corr_filter_a = xcorr(Subcarr1.*SigLOC_tot_d, PRN_advance, 'coeff');  %BOC信号与超前PRN码的互相关函数
BOCPRN_corr_filter_l = xcorr(Subcarr1.*SigLOC_tot_d, PRN_lag, 'coeff');  %BOC信号与超前PRN码的互相关函数
filter = abs(BOC_corr-0.5*(BOCPRN_corr_filter_a-BOCPRN_corr_filter_l)).^2;
%filter = Aspect + abs(BOCPRN_corr_filter_a).^2+abs(BOCPRN_corr_filter_l).^2 - abs(BOCPRN_corr_filter_a+BOCPRN_corr_filter_l).^2;
%filter = abs(BOCPRN_corr_filter_a)+abs(BOCPRN_corr_filter_l) - abs(BOCPRN_corr_filter_a+BOCPRN_corr_filter_l);
b_max = max(filter);

%导频分量
PRN_p_lag = [SigLOC_tot_p(PRN_num:end) zero_filling];        %滞后半个码片的PRN
PRN_p_advance = [a_filling SigLOC_tot_p(1:end-a_num+1)];    %超前半个码片的PRN
BOCPRN_corr_p_filter_a = xcorr(sqrt(1/11)*Subcarr2.*SigLOC_tot_p+j*sqrt(29/44)*Subcarr1.*SigLOC_tot_p, PRN_p_advance, 'coeff');  %BOC信号与PRN码的互相关函数
BOCPRN_corr_p_filter_l = xcorr(sqrt(1/11)*Subcarr2.*SigLOC_tot_p+j*sqrt(29/44)*Subcarr1.*SigLOC_tot_p, PRN_p_lag, 'coeff');  %BOC信号与超前PRN码的互相关函数
%filter_p = abs((BOC_corr_p)-0.5*(BOCPRN_corr_p_filter_a-BOCPRN_corr_p_filter_l)).^2 ;
%filter_p = abs(BOC_corr_p.*BOCPRN_corr_p_filter_a)+abs(BOC_corr_p.*BOCPRN_corr_p_filter_l);
filter_p = abs(BOCPRN_corr_p_filter_a)+abs(BOCPRN_corr_p_filter_l) - abs(BOCPRN_corr_p_filter_a+BOCPRN_corr_p_filter_l);
r_max = max(filter_p);

figure(3)
%subplot(1,2,1);
%plot(index,sqrt(BOC_corr.^2), 'b',index, filter, 'r', index, sqrt(Aspect.^2), 'm');
plot(index,sqrt(BOC_corr.^2), 'b',index, filter, 'r');
%hold on;
%plot(index,BOCPRN_corr_filter_a,'g',index, BOCPRN_corr_filter_l,'y');
xlabel('码片');
ylabel('相关函数');
axis([-1.5 1.5 -0.5 2.3]);
legend('BOC(1,1)的自相关函数','Filtered相关捕获算法');
%legend('BOC(1,1)的自相关函数','改进后的捕获算法','ASPeCT捕获算法','BOC-超前PRN互相关','BOC-滞后PRN互相关');
%title('数据分量的自相关边峰消除前后对比');
% subplot(1,2,2);
% plot(index, filter_p/r_max, 'b', index,sqrt(BOC_corr_p.^2), 'r');
% hold on;
% plot(index,BOCPRN_corr_p_filter_a,'g',index, BOCPRN_corr_p_filter_l,'k');
% xlabel('码片');
% ylabel('归一化自相关函数');
% axis([-1.5 1.5 -0.5 2.1]);
% legend('基于filter的重构相关函数','QMBOC(6,1,4/33)自相关','QMBOC-超前PRN互相关','QMBOC-滞后PRN互相关');
% title('导频分量的自相关边峰消除前后对比');

%基于filter的方法应用到BOC(6,1)信号
BOC_6_1 = SigLOC_tot_p.*Subcarr2;
BOC_6_1_corr = xcorr( BOC_6_1, 'coeff');  %BOC(6,1)信号自相关函数
BOCPRN_corr_6_1_filter_a = xcorr(BOC_6_1, PRN_p_advance, 'coeff');  %BOC信号与超前PRN码的互相关函数
BOCPRN_corr_6_1_filter_b = xcorr(BOC_6_1, PRN_p_lag, 'coeff');  %BOC信号与超前PRN码的互相关函数
%filter_6_1 = abs(BOC_6_1_corr-0.5*(BOCPRN_corr_6_1_filter_a-BOCPRN_corr_6_1_filter_b)).^2;
filter_6_1 = abs(BOCPRN_corr_6_1_filter_a).^2+abs(BOCPRN_corr_6_1_filter_b).^2 - abs(BOCPRN_corr_6_1_filter_a+BOCPRN_corr_6_1_filter_b).^2;

figure(4)
plot(index, filter_6_1, 'b', index,sqrt(Aspect_p_n.^2), 'r');
hold on;
plot(index,BOCPRN_corr_6_1_filter_a,'g',index, BOCPRN_corr_6_1_filter_b,'k');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1.5 1.5 -0.5 2.1]);
legend('基于filter的重构相关函数','ASPeCT','BOC-超前PRN互相关','BOC-滞后PRN互相关');
title('导频分量的自相关边峰消除前后对比');

%基于OQCC相关法的重构相关函数
Subcarr1_quad = sign(sin(2*pi*f_sc_a*t + pi/2));   %cos
Subcarr2_quad = sign(sin(2*pi*f_sc_b*t + pi/2));
%数据分量
QBOC_d = Subcarr1_quad.*SigLOC_tot_d;
offset_num =   1/4*1/Rc*f_sample;                                  %偏移四分之一个码片对应的时间
zero_filling_d = zeros(1, offset_num);
QBOC_d_left = [zero_filling_d QBOC_d(1:end-offset_num)];    %QBOC左偏
QBOC_d_right = [QBOC_d(offset_num+1:end) zero_filling_d];   %QBOC右偏
BOC_QBOC_d_corr_left = xcorr(Subcarr1.*SigLOC_tot_d, QBOC_d_left, 'coeff');  %BOC信号与超前QBOC信号的互相关函数
BOC_QBOC_d_corr_right = xcorr(Subcarr1.*SigLOC_tot_d, QBOC_d_right, 'coeff');  %BOC信号与滞后QBOC信号的互相关函数

OQCC_d = (abs(BOC_corr) - 0.3*abs(BOC_QBOC_d_corr_left + BOC_QBOC_d_corr_right)).^2;

%导频分量
QBOC_p = j*Subcarr1_quad.*SigLOC_tot_p + Subcarr2_quad.*SigLOC_tot_p;
offset_p_num =   1/2*1/Rc*f_sample;                                      %偏移八分之一个码片对应的时间
zero_filling_p = zeros(1, offset_p_num);
QBOC_p_left = [zero_filling_p QBOC_p(1:end-offset_p_num)];    %QBOC左偏
QBOC_p_right = [QBOC_p(offset_p_num+1:end) zero_filling_p];   %QBOC右偏
BOC_QBOC_p_corr_left = xcorr(-j*sqrt(4/33)*Subcarr2.*SigLOC_tot_p+sqrt(29/33)*Subcarr1.*SigLOC_tot_p, QBOC_p_left, 'coeff');  %BOC信号与超前QBOC信号的互相关函数
BOC_QBOC_p_corr_right = xcorr(-j*sqrt(4/33)*Subcarr2.*SigLOC_tot_p+sqrt(29/33)*Subcarr1.*SigLOC_tot_p, QBOC_p_right, 'coeff');  %BOC信号与滞后QBOC信号的互相关函数
OQCC_p = (abs(BOC_corr_p) - 0.5*abs(BOC_QBOC_p_corr_left + BOC_QBOC_p_corr_right)).^2;

figure(5)
subplot(1,2,1);
plot(index, OQCC_d, 'b', index, sqrt(Aspect.^2), 'r');
hold on;
plot(index,BOC_QBOC_d_corr_left,'g', index,BOC_QBOC_d_corr_right,'y');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1.5 1.5 -1 1.2]);
legend('基于OQCC的重构相关函数','ASPeCT','BOC-超前QBOC互相关','BOC-滞后QBOC互相关');
title('数据分量的自相关边峰消除前后对比');

subplot(1,2,2);
plot(index, OQCC_p, 'b', index,  sqrt(Aspect_p.^2), 'r');
hold on;
plot(index,BOC_QBOC_p_corr_left,'g', index,BOC_QBOC_p_corr_right,'y');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1.5 1.5 -1 1.2]);
legend('基于OQCC的重构相关函数','QMBOC(6,1,4/33)自相关','QMBOC-超前QBOC互相关','QMBOC-滞后QBOC互相关');
title('导频分量的自相关边峰消除前后对比');

%用于验证OQCC方法应用在BOC(6,1)信号中
BOC_6_1 = SigLOC_tot_p.*Subcarr2;
QBOC_6_1 = SigLOC_tot_p.*sign(sin(2*pi*f_sc_b*t + pi/2));
BOC_6_1_corr = xcorr( BOC_6_1, 'coeff');  %BOC(6,1)信号自相关函数
BOC_QBOC_6_1 = xcorr(BOC_6_1, QBOC_6_1,  'coeff');  %BOC(6,1)信号与QBOC信号互相关函数

BOC_QBOC_6_1_left = [zero_filling_p QBOC_6_1(1:end-offset_p_num)];    %QBOC左偏
BOC_QBOC_6_1_right = [QBOC_6_1(offset_p_num+1:end) zero_filling_p];   %QBOC右偏
BOC_QBOC_6_1_corr_left = xcorr(BOC_6_1, BOC_QBOC_6_1_left, 'coeff');  %BOC信号与超前QBOC信号的互相关函数
BOC_QBOC_6_1_corr_right = xcorr(BOC_6_1, BOC_QBOC_6_1_right, 'coeff');  %BOC信号与滞后QBOC信号的互相关函数

OQCC_6_1 = (abs(BOC_6_1_corr) - 0.5*abs(BOC_QBOC_6_1_corr_left + BOC_QBOC_6_1_corr_right)).^2;

figure(6)
subplot(1,2,1)
plot(index, BOC_6_1_corr, 'r', index,OQCC_6_1, 'b');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1 1 -1 1]);
legend('BOC(6,1)','BOC-QBOC');
title('BOC(6,1)信号自相关函数与BOC-QBOC信号互相关函数');
subplot(1,2,2)
plot(index, sqrt(Aspect_p_n.^2), 'r', index,OQCC_6_1, 'b');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1 1 -1 1]);
legend('BOC(6,1)','OQCC');
title('BOC(6,1)信号自相关函数与OQCC重构相关函数');

%用于观察QBOC信号偏移的距离对相关函数的影响
%offset_6_1_a =   1/8*1/Rc*f_sample;                                      %偏移八分之一个码片对应的时间
%offset_6_1_b =   floor(1/12/Rc*f_sample);
%offset_6_1_c =   1/4*1/Rc*f_sample;
%offset_6_1_d =   1/10*1/Rc*f_sample;
%zero_filling_a = zeros(1, offset_6_1_a);
%zero_filling_b = zeros(1, offset_6_1_b);
%zero_filling_c = zeros(1, offset_6_1_c);
%zero_filling_d = zeros(1, offset_6_1_d);

%BOC_QBOC_6_1_left_a = [zero_filling_a QBOC_6_1(1:end-offset_6_1_a)];    %QBOC左偏
%BOC_QBOC_6_1_right_a = [QBOC_6_1(offset_6_1_a+1:end) zero_filling_a];   %QBOC右偏
%BOC_QBOC_6_1_corr_left_a = xcorr(BOC_6_1, BOC_QBOC_6_1_left_a, 'coeff');  %BOC信号与超前QBOC信号的互相关函数
%BOC_QBOC_6_1_corr_right_a = xcorr(BOC_6_1, BOC_QBOC_6_1_right_a, 'coeff');  %BOC信号与滞后QBOC信号的互相关函数
%OQCC_6_1_a = (abs(BOC_6_1_corr) - 0.5*abs(BOC_QBOC_6_1_corr_left_a + BOC_QBOC_6_1_corr_right_a)).^2;

%BOC_QBOC_6_1_left_b = [zero_filling_b QBOC_6_1(1:end-offset_6_1_b)];    %QBOC左偏
%BOC_QBOC_6_1_right_b = [QBOC_6_1(offset_6_1_b+1:end) zero_filling_b];   %QBOC右偏
%BOC_QBOC_6_1_corr_left_b = xcorr(BOC_6_1, BOC_QBOC_6_1_left_b, 'coeff');  %BOC信号与超前QBOC信号的互相关函数
%BOC_QBOC_6_1_corr_right_b = xcorr(BOC_6_1, BOC_QBOC_6_1_right_b, 'coeff');  %BOC信号与滞后QBOC信号的互相关函数
%OQCC_6_1_b = (abs(BOC_6_1_corr) - 0.5*abs(BOC_QBOC_6_1_corr_left_b + BOC_QBOC_6_1_corr_right_b)).^2;

%BOC_QBOC_6_1_left_c = [zero_filling_c QBOC_6_1(1:end-offset_6_1_c)];    %QBOC左偏
%BOC_QBOC_6_1_right_c = [QBOC_6_1(offset_6_1_c+1:end) zero_filling_c];   %QBOC右偏
%BOC_QBOC_6_1_corr_left_c = xcorr(BOC_6_1, BOC_QBOC_6_1_left_c, 'coeff');  %BOC信号与超前QBOC信号的互相关函数
%BOC_QBOC_6_1_corr_right_c = xcorr(BOC_6_1, BOC_QBOC_6_1_right_c, 'coeff');  %BOC信号与滞后QBOC信号的互相关函数
%OQCC_6_1_c = (abs(BOC_6_1_corr) - 0.5*abs(BOC_QBOC_6_1_corr_left_c + BOC_QBOC_6_1_corr_right_c)).^2;

%BOC_QBOC_6_1_left_d = [zero_filling_d QBOC_6_1(1:end-offset_6_1_d)];    %QBOC左偏
%BOC_QBOC_6_1_right_d = [QBOC_6_1(offset_6_1_d+1:end) zero_filling_d];   %QBOC右偏
%BOC_QBOC_6_1_corr_left_d = xcorr(BOC_6_1, BOC_QBOC_6_1_left_d, 'coeff');  %BOC信号与超前QBOC信号的互相关函数
%BOC_QBOC_6_1_corr_right_d = xcorr(BOC_6_1, BOC_QBOC_6_1_right_d, 'coeff');  %BOC信号与滞后QBOC信号的互相关函数
%OQCC_6_1_d = (abs(BOC_6_1_corr) - 0.5*abs(BOC_QBOC_6_1_corr_left_d + BOC_QBOC_6_1_corr_right_d)).^2;

%figure(6)
%plot(index, OQCC_6_1_a, 'r', index,OQCC_6_1_b, 'b');
%hold on
%plot(index, OQCC_6_1_c, 'k', index,OQCC_6_1_d, 'y');
%%xlabel('码片');
%ylabel('归一化自相关函数');
%axis([-1 1 -0.2 1]);
%legend('1/8','1/2','1/4','1/10');

figure(7)
plot(index, filter_p/r_max, 'b', index,sqrt(BOC_corr_p.^2), 'r--');
xlabel('码片');
ylabel('归一化自相关函数');
axis([-1.5 1.5 -0.5 1.1]);
legend('捕获算法重构相关函数','QMBOC(6,1,4/33)自相关');
