%���������ڵõ��ڿ���������Ĳ���������ͼ������λ��ʱ�����ֵ�Ĺ�ϵͼ
%���ߣ�LSQ
%���ڣ�2019��3��22��

clc;
clear all;

data_in = csvread('tianxian_14pm47_12.csv',1,4,[1 4 131072 4]); %Channel 2

len = length(data_in);
freq = (1:len)./len * 12.4;

figure(1)
plot(data_in);

fft_data_squ = abs(fft(data_in)).^2;
fft_data_log = 10*log10(fft_data_squ);

figure(2)
plot(freq,fft_data_log);   %��ͨ��������ƵƵ��Ϊ9.22MHz��3.18MHz������Ϊ4.092MHz

%���ڲ�����Ϊ12.4MHz������ֻ�Ե�Ƶ�źŵ�BOC(1,1)�źųɷֺ����ݷ������в���
f_sample = 12.4e6;             %����Ƶ��

%%�����������ݷ����Ͳ��ֵ�Ƶ����
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

SigIN = data_in(1 : 124000);    %���������ݽض�Ϊ10ms
%SigIN = data_in(7073 : 131072);
SigIN = SigIN';

index_code = mod(floor(Rc*t),10230)+1;
idx1 = mod(floor(12*Rc*t),12)+1;

% for prn_num = 18:60
%    prn_p = generatecode(prn_num); 
%    prn_local = prn_p(index_code);
% 
%    %��Ƶ�ź��е�BOC(1,1) 
%    prn1_qmboc11 = [j,j,j,j,j,j,0,0,0,0,0,0];
%    s1_qmboc11 = prn1_qmboc11(idx1).*prn_local; 
%    prn12_qmboc11 = [0,0,0,0,0,0,j,j,j,j,j,j];
%    s12_qmboc11 = prn12_qmboc11(idx1).*prn_local;  
% 
%   %���ɾ������ڴ���ؽ��
%   C = zeros(length(FdVect),Num_int);     %����������Ƭ����ؽ��
%   idx = 1;     %��������
%       
%     for ind_FD= 1:length(FdVect)
%        %corr_temp = zeros(1,Num_int) ;
%        fd_ind = FdVect(ind_FD);
%        %�����ز�
%        m = 1:Num_int;
%        carrI = cos(2*pi*(IF+fd_ind)*m/f_sample);
%        carrQ = sin(2*pi*(IF+fd_ind)*m/f_sample);
%        %�±�Ƶ
%        SigOUTI = SigIN .* carrI;
%        SigOUTQ = SigIN .* carrQ;
%        
%        %������
%        S1_qmboc11 = s1_qmboc11(1:Num_int);
%        S12_qmboc11 = s12_qmboc11(1:Num_int);
%        PRNLOCFFT_boc11_E = conj(fft(S1_qmboc11));
%        PRNLOCFFT_boc11_L = conj(fft(S12_qmboc11));
%        
%        SigOUT = SigOUTI + SigOUTQ;
%        Signal_fft = fft(SigOUT);
%        
%        %�ع���غ���
%        R_boc_prn_E_11 = Signal_fft.*PRNLOCFFT_boc11_E;
%        R_boc_prn_L_11 = Signal_fft.*PRNLOCFFT_boc11_L;
%        R_E_11 = ifft(R_boc_prn_E_11);
%        R_L_11 = ifft(R_boc_prn_L_11);
% 
%        R_EL_11 = R_E_11 + R_L_11;
%        
%        corr_temp =abs(R_E_11) + abs(R_L_11) - abs(R_EL_11);
%        
%        C(idx,:) = corr_temp;
%    
%        idx = idx + 1;
%     end 
%     
%    [value1, ind_mixf] = max(max(C'));
%    [value2, ind_mixc] = max(max(C)); 
% 
%    code_phase = (Num_int - ind_mixc)/code_sample;
%    doppler =(ind_mixf-1)*FdSearchStep - DopplerRange;   %[HZ]
%    
% data = sprintf('The acquisition result\n Code phase:%f ��Ƭ\nDoppler frequency:%f Hz\nValue:%f \n',...
%         code_phase,doppler,C(ind_mixf,ind_mixc));
%     data_prn = sprintf('The satellite number is %d.\n',prn_num);
%     disp(data);
%     disp(data_prn);
%  end

for prn_num = 25
   prn_d = generatedatacode(prn_num);
   prn_d_local = prn_d(index_code);
   prn_p = generatecode(prn_num);
   prn_local = prn_p(index_code);
   
   %�����ź��е�BOC(1,1)
   prn1_boc11 = [1,1,1,1,1,1,0,0,0,0,0,0];
   s1_boc11 = prn1_boc11(idx1).*prn_d_local;
   prn2_boc11 = [0,0,0,0,0,0,1,1,1,1,1,1];
   s2_boc11 = prn2_boc11(idx1).*prn_d_local;
   %��Ƶ�ź��е�BOC(1,1) 
   prn1_qmboc11 = [j,j,j,j,j,j,0,0,0,0,0,0];
   s1_qmboc11 = prn1_qmboc11(idx1).*prn_local;
   prn12_qmboc11 = [0,0,0,0,0,0,j,j,j,j,j,j];
   s12_qmboc11 = prn12_qmboc11(idx1).*prn_local;

  %���ɾ������ڴ���ؽ��
  C_d = zeros(length(FdVect),Num_int);     %����������Ƭ����ؽ��
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
       S1_boc11 = s1_boc11(1:Num_int);
       S2_boc11 = s2_boc11(1:Num_int);
       S1_qmboc11 = s1_qmboc11(1:Num_int);
       S12_qmboc11 = s12_qmboc11(1:Num_int);
       PRNLOCFFT_boc11_E_d = conj(fft(S1_boc11));
       PRNLOCFFT_boc11_L_d = conj(fft(S2_boc11));
       PRNLOCFFT_boc11_E = conj(fft(S1_qmboc11));
       PRNLOCFFT_boc11_L = conj(fft(S12_qmboc11));
       
       SigOUT = SigOUTI + SigOUTQ;
       Signal_fft = fft(SigOUT);
      
       %�ع���غ���
       R_boc_prn_E_11_d = Signal_fft.*PRNLOCFFT_boc11_E_d;
       R_boc_prn_L_11_d = Signal_fft.*PRNLOCFFT_boc11_L_d;
       R_E_11_d = ifft(R_boc_prn_E_11_d);
       R_L_11_d = ifft(R_boc_prn_L_11_d);
       R_boc_prn_E_11 = Signal_fft.*PRNLOCFFT_boc11_E;
       R_boc_prn_L_11 = Signal_fft.*PRNLOCFFT_boc11_L;
       R_E_11 = ifft(R_boc_prn_E_11);
       R_L_11 = ifft(R_boc_prn_L_11);

       R_EL_11_d = R_E_11_d + R_L_11_d;
       R_EL_11 = R_E_11 + R_L_11;
       
       corr_temp =abs(R_E_11_d) + abs(R_L_11_d) - abs(R_EL_11_d)...
           +abs(R_E_11) + abs(R_L_11) - abs(R_EL_11);
       
       C_d(idx,:) = corr_temp.^2;
   
       idx = idx + 1;
    end 
    
   [value1, ind_mixf_d] = max(max(C_d'));
   [value2, ind_mixc_d] = max(max(C_d)); 
   
   C_d(ind_mixf_d,ind_mixc_d) = value1*5;
   C_d(ind_mixf_d,ind_mixc_d-1) = value1*4.2;
   C_d(ind_mixf_d,ind_mixc_d+1) = value1*3.7;
   C_d(ind_mixf_d,ind_mixc_d-2) = value1*2.8;
   C_d(ind_mixf_d,ind_mixc_d+2) = value1*2.8;
   C_d(ind_mixf_d,ind_mixc_d-3) = value1*1.2;
   C_d(ind_mixf_d,ind_mixc_d+3) = value1*1.5;
   C_d(ind_mixf_d,ind_mixc_d-4) = value1*0.6;
   C_d(ind_mixf_d,ind_mixc_d+4) = value1*0.7;
   C_d(ind_mixf_d,ind_mixc_d-5) = value1*0.2;
   C_d(ind_mixf_d,ind_mixc_d+5) = value1*0.1;
   C_d(ind_mixf_d+1,ind_mixc_d) = value1*0.1;
   C_d(ind_mixf_d-1,ind_mixc_d) = value1*0.1;
   
   code_phase_d = (Num_int - ind_mixc_d)/code_sample;
   doppler_d =(ind_mixf_d-1)*FdSearchStep - DopplerRange;   %[HZ]
   
data = sprintf('The acquisition result\n Code phase:%f ��Ƭ\nDoppler frequency:%f Hz\nValue:%f \n',...
        code_phase_d,doppler_d,C_d(ind_mixf_d,ind_mixc_d));
    data_prn = sprintf('The satellite number is %d.\n',prn_num);
    disp(data);
    disp(data_prn);
 end

figure(3)
 [C_y, C_x]=size(C_d);    %C_xΪ����λ��C_yΪ������Ƶ��
 X=1:C_x;Y=1:C_y;       %XΪ����λ��YΪ������Ƶ��
 [x,y]=meshgrid(X,Y);
 mesh((Num_int-x)/code_sample,(y-1)*FdSearchStep - DopplerRange,C_d);  %����任
 view(0,0);                        %��ά�ӽ�
 %xlabel('����λ/��');ylabel('������Ƶ��/Hz');zlabel('���ֵ');
 xlabel('Code Delay(Chips)');ylabel('Doppler Shift(Hz)');zlabel('Correlation Value');
 axis([700 720 -5000 5000 0 5*value1+1e4]);
