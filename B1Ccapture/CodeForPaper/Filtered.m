%Author:LSQ
%Date:2019/4
%Description: 论文3.2.6章节.

clc;
close all;

set(0,'defaultfigurecolor','w'); %将仿真图背景设置为白色

%验证Filtered捕获算法的正确性
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

num_sample = floor(f_sample/Rc);

[BOCPRN_corr lag] = xcorr(Subcarr1.*SigLOC_tot_d, SigLOC_tot_d, 'coeff');  %BOC信号与PRN码的互相关函数
BOC_corr = xcorr(Subcarr1.*SigLOC_tot_d, 'coeff');    %BOC的自相关函数
%横坐标尺度变换
index = lag / num_sample;

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
b_max = max(filter);

%导频分量
PRN_p_lag = [SigLOC_tot_p(PRN_num:end) zero_filling];        %滞后半个码片的PRN
PRN_p_advance = [a_filling SigLOC_tot_p(1:end-a_num+1)];    %超前半个码片的PRN
BOCPRN_corr_p_filter_a = xcorr(sqrt(1/11)*Subcarr2.*SigLOC_tot_p+j*sqrt(29/44)*Subcarr1.*SigLOC_tot_p, PRN_p_advance, 'coeff');  %BOC信号与PRN码的互相关函数
BOCPRN_corr_p_filter_l = xcorr(sqrt(1/11)*Subcarr2.*SigLOC_tot_p+j*sqrt(29/44)*Subcarr1.*SigLOC_tot_p, PRN_p_lag, 'coeff');  %BOC信号与超前PRN码的互相关函数
filter_p = abs(BOCPRN_corr_p_filter_a)+abs(BOCPRN_corr_p_filter_l) - abs(BOCPRN_corr_p_filter_a+BOCPRN_corr_p_filter_l);
r_max = max(filter_p);

figure(1)
plot(index,sqrt(BOC_corr.^2), 'b',index, filter, 'r');
xlabel('相位延时(码片)');
ylabel('相关值');
axis([-1.5 1.5 -0.5 2.3]);
legend('BOC(1,1)的自相关函数','Filtered相关捕获算法');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

