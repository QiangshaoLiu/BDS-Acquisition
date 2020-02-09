%Author:LSQ
%Date:2019/3
%Description:This program is for Beidou B1C satellite signal acquisition.
%Parameters: 
%          Sampling rate: 12.4MHz
%          Channel 1: Beidou B3
%          Channel 2: GPS L1
%          Sampling depth: 131072
%     ASPeCT算法
clc;
close all;

data_in_1 = csvread('tianxian_14pm47_12.csv',1,3,[1 3 131072 3]);%Channel 1
data_in_2 = csvread('tianxian_14pm47_12.csv',1,4,[1 4 131072 4]);%Channel 2

len = length(data_in_2);
freq = (1:len)./len * 12.4;

figure(1)
plot(data_in_2);

fft_data_squ = abs(fft(data_in_2)).^2;
fft_data_log = 10*log10(fft_data_squ);

figure(2)
plot(freq,fft_data_log);   %带通采样，中频频点为9.22MHz和3.18MHz，带宽为4.092MHz

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
FdSearchStep = 40;      %[Hz]
DopplerRange = 5000;      %[Hz]
code_sample = floor(f_sample/Rc);   %单个码片所对应的采样数
FdVect= -DopplerRange:FdSearchStep:DopplerRange;     %多普勒频移搜索范围

Subcarr1 = sign(sin(2*pi*f_sc_a*t));

SigIN = data_in_2(1 : 124000);    %将卫星数据截短为10ms
SigIN = SigIN';

for prn_num = 23
   prn_p = generatecode(prn_num);
   index_code = mod(floor(Rc*t),10230)+1;
   prn_local = prn_p(index_code);

%导频信号中的BOC(1,1)
   B1C_poilt = j*prn_local.*Subcarr1;
   BOCLOCFFT_boc11 = conj(fft(B1C_poilt));

  %生成矩阵用于存相关结果
  C = zeros(length(FdVect),Num_int);     %用于所有码片的相关结果
  idx = 1;     %矩阵行数
  
    for ind_FD= 1:length(FdVect)
       %corr_temp = zeros(1,Num_int) ;
       fd_ind = FdVect(ind_FD);
       %本地载波
       m = 1:Num_int;
       carrI = cos(2*pi*(IF+fd_ind)*m/f_sample);
       carrQ = sin(2*pi*(IF+fd_ind)*m/f_sample);
       %下变频
       SigOUTI = SigIN .* carrI;
       SigOUTQ = SigIN .* carrQ;
       
       %本地码
       PRNLOCFFT_boc11 = conj(fft(prn_local));
       
       SigOUT = SigOUTI + SigOUTQ;
       Signal_fft = fft(SigOUT);
       
       %重构相关函数
       R_boc_prn =ifft( Signal_fft.*PRNLOCFFT_boc11 );
       R_boc_boc =ifft( Signal_fft.*BOCLOCFFT_boc11 );
       
       corr_temp =abs(R_boc_boc).^2 - abs(R_boc_prn).^2;
       
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

figure(3)
[C_y, C_x]=size(C);    %C_x为码相位，C_y为多普勒频移
X=1:C_x;Y=1:C_y;       %X为码相位，Y为多普勒频移
[x,y]=meshgrid(X,Y);
mesh((Num_int-x)/code_sample,(y-1)*FdSearchStep - DopplerRange,C);  %坐标变换
view(0,0);                        %二维视角
xlabel('码相位/码');ylabel('多普勒频移/Hz');zlabel('相关值');
axis([4690 4710 -3700 -2700 0 value1+1e9]);

