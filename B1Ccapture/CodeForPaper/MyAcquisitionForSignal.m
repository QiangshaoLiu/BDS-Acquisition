%Author:LSQ
%Date:2019/4
%Description: 利用学位论文设计的B1C信号捕获算法对采集到的卫星信号数据进行处理，得到接收机的收星情况.
%Parameters: 
%          Sampling rate: 12.4MHz
%          USB3.0 AD 2: GPS L1
%          Sampling depth: 131072

clc;
close all;

set(0,'defaultfigurecolor','w'); %将仿真图背景设置为白色

[data_num data_order data_AD_in] = textread('201904031204.txt','%d%d%d');


%由于采样率为12.4MHz，所以仅对导频信号的BOC(1,1)信号成分进行捕获
f_sample = 12.4e6;             %采样频率

%%产生本地导频信号
%基本参数
f_sc_a = 1.023e6 ;                 %BOC(1,1)子载波速率
f_sc_b = 6*1.023e6 ;               %BOC(6,1)子载波速率
Rc = 1.023e6;                      %主码码速率
T_process = 10e-3;                 %处理时间
T_int = 10e-3;                      %相关运算时间
t = 0 : 1/f_sample : T_process - 1/f_sample;
n = 0:length(t)-1;                 
j=sqrt(-1);
pi = 3.141592654;                  %圆周率
Num_int = floor(f_sample * T_int); %相干积分时间所对应的采样点数
IF = 3.18e6;            %[Hz]
FdSearchStep = 50;      %[Hz]
DopplerRange = 5000;      %[Hz]
code_sample = floor(f_sample/Rc);   %单个码片所对应的采样数
FdVect= -DopplerRange:FdSearchStep:DopplerRange;     %多普勒频移搜索范围

SigIN = data_AD_in(2 : 124001);    %将卫星数据截短为10ms
SigIN = SigIN';

for prn_num = 18:60
   prn_p = generatecode(prn_num);
   prn_d = generatedatacode(prn_num);
   
   index_code = mod(floor(Rc*t),10230)+1;
   prn_local = prn_p(index_code);
   prn_d_local = prn_d(index_code);

  %导频信号中的BOC(1,1)
   idx1 = mod(floor(12*Rc*t),12)+1;
   prn1_qmboc11 = [j,j,j,j,j,j,0,0,0,0,0,0];
   s1_qmboc11 = prn1_qmboc11(idx1).*prn_local;
   
   prn12_qmboc11 = [0,0,0,0,0,0,j,j,j,j,j,j];
   s12_qmboc11 = prn12_qmboc11(idx1).*prn_local;

   %数据信号中的BOC(1,1)
   prn1_boc11 = [1,1,1,1,1,1,0,0,0,0,0,0];
   s1_boc11 = prn1_boc11(idx1).*prn_d_local;
   
   prn2_boc11 = [0,0,0,0,0,0,1,1,1,1,1,1];
   s2_boc11 = prn2_boc11(idx1).*prn_d_local;
   
  %生成矩阵用于存相关结果
  C = zeros(length(FdVect),Num_int);     %用于所有码片的相关结果
  idx = 1;     %矩阵行数
      
    for ind_FD= 1:length(FdVect)
       fd_ind = FdVect(ind_FD);
       %本地载波
       m = 1:Num_int;
       carrI = cos(2*pi*(IF+fd_ind)*m/f_sample);
       carrQ = sin(2*pi*(IF+fd_ind)*m/f_sample);
       %下变频
       SigOUTI = SigIN .* carrI;
       SigOUTQ = SigIN .* carrQ;
       
       %本地码
       S1_qmboc11 = s1_qmboc11(1:Num_int);
       S12_qmboc11 = s12_qmboc11(1:Num_int);
       S1_boc11 = s1_boc11(1:Num_int);
       S2_boc11 = s2_boc11(1:Num_int);
       
       PRNLOCFFT_boc11_E = conj(fft(S1_qmboc11));
       PRNLOCFFT_boc11_L = conj(fft(S12_qmboc11));
       PRNLOCFFT_boc11_d_E = conj(fft(S1_boc11));
       PRNLOCFFT_boc11_d_L = conj(fft(S2_boc11));
       
       SigOUT = SigOUTI + SigOUTQ;
       Signal_fft = fft(SigOUT);
       
       %重构相关函数
       R_boc_prn_E_11 = Signal_fft.*PRNLOCFFT_boc11_E;
       R_boc_prn_L_11 = Signal_fft.*PRNLOCFFT_boc11_L;
       R_boc_prn_d_E_11 = Signal_fft.*PRNLOCFFT_boc11_d_E;
       R_boc_prn_d_L_11 = Signal_fft.*PRNLOCFFT_boc11_d_L;
       R_E_11 = ifft(R_boc_prn_E_11);
       R_L_11 = ifft(R_boc_prn_L_11);
       R_d_E_11 = ifft(R_boc_prn_d_E_11);
       R_d_L_11 = ifft(R_boc_prn_d_L_11);

       R_EL_11 = R_E_11 + R_L_11;
       R_d_EL_11 = R_d_E_11 + R_d_L_11;
       
       corr_temp =abs(R_E_11) + abs(R_L_11) - abs(R_EL_11)...
           + abs(R_d_E_11) + abs(R_d_L_11) - abs(R_d_EL_11);
       
       C(idx,:) = corr_temp;
   
       idx = idx + 1;
    end 
    
   [value1, ind_mixf] = max(max(C'));
   [value2, ind_mixc] = max(max(C)); 

   code_phase = (Num_int - ind_mixc)/code_sample;
   doppler =(ind_mixf-1)*FdSearchStep - DopplerRange;   %[HZ]
   
data = sprintf('The acquisition result\n Code phase:%f 码片\nDoppler frequency:%f Hz\nValue:%f \n',...
        code_phase,doppler,C(ind_mixf,ind_mixc));
    data_prn = sprintf('The satellite number is %d.\n',prn_num);
    disp(data);
    disp(data_prn);
  
end

 
