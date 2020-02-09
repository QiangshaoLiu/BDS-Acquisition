%本程序给出了使用伪相关函数法得到的相关函数波形图
%作者：LSQ
%日期：2019年3月21日
clc;
close all;

f_sample = 30*1.023e6;             %采样频率
f_sc_a = 1.023e6 ;                 %数据分量子载波速率
f_sc_b = 6*1.023e6 ;               %导频分量子载波速率
Rc = 1.023e6;                      %主码码速率
T_process = 10e-3;                 %处理时间
t = 0 : 1/f_sample : T_process - 1/f_sample;                
j=sqrt(-1);
pi = 3.141592654;                  %圆周率

%%模拟产生接收信号
subcarr1 = sign(sin(2*pi*f_sc_a*t));
subcarr1(1) = 1;
subcarr2 = sign(sin(2*pi*f_sc_b*t));
subcarr2(1) = 1;
code_r = generatecode(1);           %接收信号由PRN=1的扩频码序列调制
codeSample_r = code_r(mod(floor(t*Rc),10230)+1);
Qmboc_p = sqrt(1/11)*codeSample_r.*subcarr2 + ...
    j*sqrt(29/44)*codeSample_r.*subcarr1;

BOC_6_1 = codeSample_r.*subcarr2;
BOC_1_1 = codeSample_r.*subcarr1;

%产生本地测距码序列
prn_p = generatecode(1);
index_code = mod(floor(Rc*t),10230)+1;
prn_local = prn_p(index_code);

%%导频信号QMBOC(6,1,4/33)
idx1 = mod(floor(12*Rc*t),12)+1;
prn1_qmboc11 = [j*sqrt(29/44),0,0,0,0,0,0,0,0,0,0,0];
s1_qmboc11 = prn1_qmboc11(idx1).*prn_local;
[g1_qmboc11 x]= xcorr(Qmboc_p, s1_qmboc11, 'coeff');
prn1_qmboc61 = [sqrt(1/11),sqrt(1/11),sqrt(1/11),sqrt(1/11),...
    sqrt(1/11),sqrt(1/11),0,0,0,0,0,0];
s1_qmboc61 = prn1_qmboc61(idx1).*prn_local;
g1_qmboc61 = xcorr(Qmboc_p, s1_qmboc61, 'coeff');
prn12_qmboc11 = [0,0,0,0,0,0,0,0,0,0,0,j*sqrt(29/44)];
s12_qmboc11 = prn12_qmboc11(idx1).*prn_local;
g12_qmboc11 = xcorr(Qmboc_p, s12_qmboc11, 'coeff');
prn12_qmboc61 = [0,0,0,0,0,0,sqrt(1/11),sqrt(1/11),sqrt(1/11),sqrt(1/11),sqrt(1/11),sqrt(1/11)];
s12_qmboc61 = prn12_qmboc61(idx1).*prn_local;
g12_qmboc61 = xcorr(Qmboc_p, s12_qmboc61, 'coeff');

corr_sum_qmboc = abs(g1_qmboc11)+abs(g1_qmboc61)+abs(g12_qmboc11)+abs(g12_qmboc61)...
    -abs(g1_qmboc11+g12_qmboc11)-abs(g1_qmboc61+g12_qmboc61);
max_num = max(corr_sum_qmboc);
corr_qmboc = xcorr(Qmboc_p, Qmboc_p, 'coeff'); 

%%BOC(1,1)信号自相关边峰消除
prn1_boc11 = [1,1,1,1,1,1,0,0,0,0,0,0];
s1_boc11 = prn1_boc11(idx1).*prn_local;
g1_boc11 = xcorr(BOC_1_1, s1_boc11, 'coeff');

prn2_boc11 = [0,0,0,0,0,0,1,1,1,1,1,1];
s2_boc11 = prn2_boc11(idx1).*prn_local;
g2_boc11 = xcorr(BOC_1_1, s2_boc11, 'coeff');

corr_sum_boc11 = abs(g1_boc11)+abs(g2_boc11)-abs(g1_boc11+g2_boc11);
corr_boc11 = xcorr(BOC_1_1, BOC_1_1, 'coeff'); 

%%BOC(6,1)信号自相关边峰消除
prn1_boc61 = [1,0,0,0,0,0,0,0,0,0,0,0];
s1_boc61 = prn1_boc61(idx1).*prn_local;
g1_boc61 = xcorr(BOC_6_1, s1_boc61, 'coeff');
prn2_boc61 = [0,0,0,0,0,0,0,0,0,0,0,1];
s2_boc61 = prn2_boc61(idx1).*prn_local;
g2_boc61 = xcorr(BOC_6_1, s2_boc61, 'coeff');
prn3_boc61 = [0,1,0,0,0,0,0,0,0,0,0,0];
s3_boc61 = prn3_boc61(idx1).*prn_local;
g3_boc61 = xcorr(BOC_6_1, s3_boc61, 'coeff');
prn4_boc61 = [0,0,0,0,0,0,0,0,0,0,1,0];
s4_boc61 = prn4_boc61(idx1).*prn_local;
g4_boc61 = xcorr(BOC_6_1, s4_boc61, 'coeff');
prn5_boc61 = [0,0,1,0,0,0,0,0,0,0,0,0];
s5_boc61 = prn5_boc61(idx1).*prn_local;
g5_boc61 = xcorr(BOC_6_1, s5_boc61, 'coeff');
prn6_boc61 = [0,0,0,0,0,0,0,0,0,1,0,0];
s6_boc61 = prn6_boc61(idx1).*prn_local;
g6_boc61 = xcorr(BOC_6_1, s6_boc61, 'coeff');
prn7_boc61 = [0,0,0,1,0,0,0,0,0,0,0,0];
s7_boc61 = prn7_boc61(idx1).*prn_local;
g7_boc61 = xcorr(BOC_6_1, s7_boc61, 'coeff');
prn8_boc61 = [0,0,0,0,0,0,0,0,1,0,0,0];
s8_boc61 = prn8_boc61(idx1).*prn_local;
g8_boc61 = xcorr(BOC_6_1, s8_boc61, 'coeff');
prn9_boc61 = [0,0,0,0,1,0,0,0,0,0,0,0];
s9_boc61 = prn9_boc61(idx1).*prn_local;
g9_boc61 = xcorr(BOC_6_1, s9_boc61, 'coeff');
prn10_boc61 = [0,0,0,0,0,0,0,1,0,0,0,0];
s10_boc61 = prn10_boc61(idx1).*prn_local;
g10_boc61 = xcorr(BOC_6_1, s10_boc61, 'coeff');
prn11_boc61 = [0,0,0,0,0,1,0,0,0,0,0,0];
s11_boc61 = prn11_boc61(idx1).*prn_local;
g11_boc61 = xcorr(BOC_6_1, s11_boc61, 'coeff');
prn12_boc61 = [0,0,0,0,0,0,1,0,0,0,0,0];
s12_boc61 = prn12_boc61(idx1).*prn_local;
g12_boc61 = xcorr(BOC_6_1, s12_boc61, 'coeff');

corr_sum_boc61 = abs(g1_boc61)+abs(g2_boc61) - abs(g1_boc61+g2_boc61);
corr_boc61 = xcorr(BOC_6_1, BOC_6_1, 'coeff'); 

%横坐标尺度变换
index = x / floor(f_sample/Rc);

 %导频分量自相关函数前后对比
 figure(1)
 plot(index, corr_qmboc,'b',index,corr_sum_qmboc,'m');
 xlabel('码片');
 ylabel('归一化自相关函数');
 axis([-1.5 1.5 -0.5 2.1]);
 %数据分量自相关函数前后对比
 figure(2)
 plot(index, corr_sum_boc11,'b',index,abs(corr_boc11),'m');
 %legend('BOC(1,1)自相关','伪相关函数');
 %xlabel('码片');
 %ylabel('归一化自相关函数');
 legend('BOC(1,1) Autocorrelation','Pseudo-Correlation Function');
 xlabel('Code Delay(Chips)');
 ylabel('Correlation Function');
 axis([-1.5 1.5 -0.2 1.5]);
 grid on;
 %导频分量自相关函数前后对比和BOC(6,1)伪相关函数
 figure(3)
 subplot(2,1,1);
 plot(index,corr_boc61,'b',index, corr_sum_boc61,'m');
 legend('BOC(6,1)自相关','伪相关函数');
 xlabel('码片');
 ylabel('归一化自相关函数');
 axis([-1.5 1.5 -1 1]);
 subplot(2,1,2);
 plot(index, corr_qmboc,'b',index,corr_sum_qmboc,'m');
 legend('QMBOC(6,1,4/33)自相关','伪相关函数');
 xlabel('码片');
 ylabel('归一化自相关函数');
 axis([-1.5 1.5 -0.5 2.1]);

 bocprncorr = xcorr(BOC_1_1,prn_local,'coeff');
 figure(4)
 plot(index, g1_boc11,'b',index,g2_boc11,'r');
 legend('R1','R2');
 %xlabel('码片');
 %ylabel('归一化自相关函数');
 xlabel('Code Delay(Chips)');
 ylabel('Correlation Function');
 axis([-1.5 1.5 -1 1]);
 grid on;
 
 figure(5)
 plot(index,corr_boc61,'b',index, corr_sum_boc61,'m');
 legend('BOC(6,1)Autocorrelation','Pseudo-Correlation Function');
 xlabel('Code Delay(Chips)');
 ylabel('Correlation Function');
 axis([-1.5 1.5 -0.8 1.2]);
 grid on;
 
 figure(6)
 plot(index, corr_qmboc,'b',index,corr_sum_qmboc,'m');
 legend('QMBOC(6,1,4/33)Autocorrelation','Pseudo-Correlation Function');
 xlabel('Code Delay(Chips)');
 ylabel('Correlation Function');
 axis([-1.5 1.5 -0.5 2.1]);
 grid on;