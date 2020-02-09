%本程序用于验证PCF算法
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

f_sample = 36*1.023e6;             %采样频率
f_sc_a = 1.023e6 ;                 %数据分量子载波速率
f_sc_b = 6*1.023e6 ;               %导频分量子载波速率
T_process = 10e-3;                 %处理时间

t = 0 : 1/f_sample : T_process - 1/f_sample;

j=sqrt(-1);

num = f_sample*T_process;   %总采样点数
Subcarr1 = zeros(1,num);
Subcarr2 = zeros(1,num);

Rc = 1.023e6;       %主码码速率
Tc = f_sample/Rc;
M = 2;       %调制阶数
Ts = Tc/M;   %分段

for i = 1 : num
    if mod(i,Tc)<Ts
        Subcarr1(i)=1;
    else
        Subcarr1(i)=-1;
    end
end

for i = 1 : num
    if mod(i,Tc/6) < Tc/12
        Subcarr2(i)=1;
    else
        Subcarr2(i)=-1;
    end
end

n=0:length(t)-1;
ind_cod = mod(floor(n*Rc/f_sample),WeilCodelength)+1;
SigLOC_tot_d = WeilCode(ind_cod);
SigLOC_tot_p = WeilCode_p(ind_cod);

%BOC(1,1)信号
BOC_1_1 = Subcarr1.*SigLOC_tot_p;

S1 = zeros(1,num);

for i = 1 : num
    if mod(i,Tc)<Ts
        S1(i)=1;
    else
        S1(i)=0;
    end
end

S2 = zeros(1,num);

for i = 1 : num
    if mod(i,Tc)<Ts
        S2(i)=0;
    else
        S2(i)=-1;
    end
end

Signal_1 = S1.*SigLOC_tot_p;
Signal_2 = S2.*SigLOC_tot_p;

[R_s1_boc lag] = xcorr(BOC_1_1, Signal_1, 'coeff');
R_s2_boc = xcorr(BOC_1_1, Signal_2, 'coeff');
A_boc = xcorr(BOC_1_1, 'coeff');
A_boc_prn = xcorr(BOC_1_1, SigLOC_tot_p,'coeff');

R = A_boc.^2 - A_boc_prn.^2;
R1 = abs(R_s1_boc) + abs(R_s2_boc) - abs(R_s1_boc + R_s2_boc);

num_sample = floor(f_sample/Rc);

%横坐标尺度变换
index = lag / num_sample;
figure(1)
plot(index, R,'b', index, R1,'r');
axis([-1.5 1.5 -0.5 1.5]);

%BOC(6,1)信号
BOC_6_1 = Subcarr2.*SigLOC_tot_p;

S3 = zeros(1,num);

for i = 1 : num
    if mod(i,Tc) < Tc/12
        S3(i)=11;
    elseif mod(i,Tc) >= Tc*11/12
        S3(i)=1;
    else
        S3(i)=0;
    end
end

S4 = zeros(1,num);

for i = 1 : num
    if mod(i,Tc) < Tc/12
        S4(i)=1;
    elseif   mod(i,Tc) >= 11*Tc/12 
        S4(i)=11;
    else
        S4(i)=0;    
    end
end

Signal_3 = S3.*SigLOC_tot_p;
Signal_4 = S4.*SigLOC_tot_p;

R_s3_boc = xcorr(BOC_6_1, Signal_3, 'coeff');
R_s4_boc = xcorr(BOC_6_1, Signal_4, 'coeff');
A_boc_1 = xcorr(BOC_6_1, 'coeff');
A_boc_prn_1 = xcorr(BOC_6_1, SigLOC_tot_p,'coeff');

R2 = A_boc_1.^2 - A_boc_prn_1.^2;
R3 = abs(R_s3_boc) + abs(R_s4_boc) - abs(R_s3_boc - R_s4_boc) ;
figure(2)
plot(index, sqrt(R2).^2,'b', index, R3,'r');
axis([-1.5 1.5 -0.5 1.5]);


