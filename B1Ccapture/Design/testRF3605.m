%�ó������ڲ���3605�޸����ú��ץ��Ч����3605�������£�
%ͨ��һ��ΪL1Ƶ�㣬�˲�������Ϊ20.46MHz,��ƵƵ��Ϊ15.58MHz��
%ͨ����ΪB1Ƶ�㣬�˲�������Ϊ4.092MHz����ƵƵ��Ϊ15.902MHz��
%���Ϊģ������������������Ȼ������

%%
%clc
%close all
id0 = 1;
id1 = 1;
%data_in_L1 = sig_AD1_IN90(id0:end,1);
%data_in_B1 = sig_AD2_IN90(id1:end,1);
%data_in_L1 = AD1_IN_IBUF90(id0:end,1);
data_in_B1 = AD2_IN_IBUF90(id1:end,1);

len = length(data_in_B1);

fs = 12.4;

freq = (1:len)./len * fs;   %��MHzΪ��λ

%fft_data_squ_L1 = abs(fft(data_in_L1)).^2;    %����
fft_data_squ_B1 = abs(fft(data_in_B1)).^2; 

%fft_data_log_L1 = 10*log10(fft_data_squ_L1);
fft_data_log_B1 = 10*log10(fft_data_squ_B1);

figure(1);
subplot(2,1,1);
%plot(data_in_L1);
subplot(2,1,2);
plot(data_in_B1);

figure(2);
subplot(2,1,1);
%plot(freq,fft_data_log_L1);
xlabel('MHz');ylabel('dB');
subplot(2,1,2);
plot(freq,fft_data_log_B1);
xlabel('MHz');ylabel('dB');

%FreqId_L1_max = min(find( fft_data_squ_L1( 1:floor(len/2) )== max( fft_data_squ_L1( 1:floor(len/2) ) ) ));
%FreqId_B1_max = min(find( fft_data_squ_B1(1:floor(len/2))== max(fft_data_squ_B1(1:floor(len/2))) ));
%FreqId_L1_max = find( fft_data_squ_L1( 1:floor(len/2) )== max( fft_data_squ_L1( 1:floor(len/2) ) ) );
FreqId_B1_max = find( fft_data_squ_B1(1:floor(len/2))== max(fft_data_squ_B1(1:floor(len/2))) );

%FreqId0_L1 = round( (15.58/5-10.23)/fs*len );
%FreqId1_L1 = round( (15.58/5+10.23)/fs*len );
FreqId0_B1 = round( (15.902/5-2.046)/fs*len );
FreqId1_B1 = round( (15.902/5+2.046)/fs*len );

f_error = round(0.004/fs*len);   %Ƶƫ���

%Noi_L1 = sum(fft_data_squ_L1(FreqId0_L1:FreqId_L1_max-f_error-1)) + sum(fft_data_squ_L1(FreqId_L1_max+f_error+1:FreqId1_L1));
Noi_B1 = sum(fft_data_squ_B1(FreqId0_B1:FreqId_B1_max-f_error-1)) + sum(fft_data_squ_B1(FreqId_B1_max+f_error+1:FreqId1_B1));

%S_L1 = sum(fft_data_squ_L1(FreqId_L1_max-f_error:FreqId_L1_max+f_error));
S_B1 = sum(fft_data_squ_B1(FreqId_B1_max-f_error:FreqId_B1_max+f_error));

%Snr_C1 = 10*log10(S_L1/Noi_L1);
Snr_C2 = 10*log10(S_B1/Noi_B1);

%f1 = FreqId_L1_max/len * fs;
f2 = FreqId_B1_max/len * fs;

disp('ͨ��һ����Ƶ�ʣ�');
%disp(['f1 = ',num2str(f1)]);

disp('ͨ��������Ƶ�ʣ�');
disp(['f2 = ',num2str(f2)]);

disp('ͨ��һ����ȣ�');
%disp(['Snr_C1 = ',num2str(Snr_C1)]);

disp('ͨ��������ȣ�');
disp(['Snr_C2 = ',num2str(Snr_C2)]);


