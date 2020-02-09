%Author:LSQ
%Date:2019/3
%Description:This program is for Beidou B1C satellite signal acquisition.
%Parameters: 
%          Sampling rate: 12.4MHz
%          Channel 1: Beidou B3
%          Channel 2: GPS L1
%          Sampling depth: 131072
%     ASPeCT�㷨
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
plot(freq,fft_data_log);   %��ͨ��������ƵƵ��Ϊ9.22MHz��3.18MHz������Ϊ4.092MHz

%���ڲ�����Ϊ12.4MHz�����Խ��Ե�Ƶ�źŵ�BOC(1,1)�źųɷֽ��в���
f_sample = 12.4e6;             %����Ƶ��

%%�������ص�Ƶ�ź�
%��������
f_sc_a = 1.023e6 ;                 %BOC(1,1)���ز�����
f_sc_b = 6*1.023e6 ;               %BOC(6,1)���ز�����
Rc = 1.023e6;                      %����������
T_process = 10e-3;                 %����ʱ��
T_int = 10e-3;                      %�������ʱ��
t = 0 : 1/f_sample : T_process - 1/f_sample;
n = 0:length(t)-1;                 
j=sqrt(-1);
pi = 3.141592654;                  %Բ����
Num_int = floor(f_sample * T_int); %��ɻ���ʱ������Ӧ�Ĳ�������
IF = 3.18e6;            %[Hz]
FdSearchStep = 40;      %[Hz]
DopplerRange = 5000;      %[Hz]
code_sample = floor(f_sample/Rc);   %������Ƭ����Ӧ�Ĳ�����
FdVect= -DopplerRange:FdSearchStep:DopplerRange;     %������Ƶ��������Χ

Subcarr1 = sign(sin(2*pi*f_sc_a*t));

SigIN = data_in_2(1 : 124000);    %���������ݽض�Ϊ10ms
SigIN = SigIN';

for prn_num = 23
   prn_p = generatecode(prn_num);
   index_code = mod(floor(Rc*t),10230)+1;
   prn_local = prn_p(index_code);

%��Ƶ�ź��е�BOC(1,1)
   B1C_poilt = j*prn_local.*Subcarr1;
   BOCLOCFFT_boc11 = conj(fft(B1C_poilt));

  %���ɾ������ڴ���ؽ��
  C = zeros(length(FdVect),Num_int);     %����������Ƭ����ؽ��
  idx = 1;     %��������
  
    for ind_FD= 1:length(FdVect)
       %corr_temp = zeros(1,Num_int) ;
       fd_ind = FdVect(ind_FD);
       %�����ز�
       m = 1:Num_int;
       carrI = cos(2*pi*(IF+fd_ind)*m/f_sample);
       carrQ = sin(2*pi*(IF+fd_ind)*m/f_sample);
       %�±�Ƶ
       SigOUTI = SigIN .* carrI;
       SigOUTQ = SigIN .* carrQ;
       
       %������
       PRNLOCFFT_boc11 = conj(fft(prn_local));
       
       SigOUT = SigOUTI + SigOUTQ;
       Signal_fft = fft(SigOUT);
       
       %�ع���غ���
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
   
data = sprintf('The acquisition result\n Code phase:%f ��Ƭ\nDoppler frequency:%f Hz\nValue:%f \n',...
        code_phase,doppler,C(ind_mixf,ind_mixc));
    data_prn = sprintf('The satellite number is %d.\n',prn_num);
    disp(data);
    disp(data_prn);
  
 end

figure(3)
[C_y, C_x]=size(C);    %C_xΪ����λ��C_yΪ������Ƶ��
X=1:C_x;Y=1:C_y;       %XΪ����λ��YΪ������Ƶ��
[x,y]=meshgrid(X,Y);
mesh((Num_int-x)/code_sample,(y-1)*FdSearchStep - DopplerRange,C);  %����任
view(0,0);                        %��ά�ӽ�
xlabel('����λ/��');ylabel('������Ƶ��/Hz');zlabel('���ֵ');
axis([4690 4710 -3700 -2700 0 value1+1e9]);

