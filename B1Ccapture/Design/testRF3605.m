%该程序用于测试3605修改配置后的抓数效果，3605配置如下：
%通道一改为L1频点，滤波器带宽为20.46MHz,中频频率为15.58MHz；
%通道二为B1频点，滤波器带宽为4.092MHz，中频频率为15.902MHz；
%输出为模拟差分输出，数字输出仍然保留。

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

freq = (1:len)./len * fs;   %以MHz为单位

%fft_data_squ_L1 = abs(fft(data_in_L1)).^2;    %求功率
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

f_error = round(0.004/fs*len);   %频偏误差

%Noi_L1 = sum(fft_data_squ_L1(FreqId0_L1:FreqId_L1_max-f_error-1)) + sum(fft_data_squ_L1(FreqId_L1_max+f_error+1:FreqId1_L1));
Noi_B1 = sum(fft_data_squ_B1(FreqId0_B1:FreqId_B1_max-f_error-1)) + sum(fft_data_squ_B1(FreqId_B1_max+f_error+1:FreqId1_B1));

%S_L1 = sum(fft_data_squ_L1(FreqId_L1_max-f_error:FreqId_L1_max+f_error));
S_B1 = sum(fft_data_squ_B1(FreqId_B1_max-f_error:FreqId_B1_max+f_error));

%Snr_C1 = 10*log10(S_L1/Noi_L1);
Snr_C2 = 10*log10(S_B1/Noi_B1);

%f1 = FreqId_L1_max/len * fs;
f2 = FreqId_B1_max/len * fs;

disp('通道一中心频率：');
%disp(['f1 = ',num2str(f1)]);

disp('通道二中心频率：');
disp(['f2 = ',num2str(f2)]);

disp('通道一信噪比：');
%disp(['Snr_C1 = ',num2str(Snr_C1)]);

disp('通道二信噪比：');
disp(['Snr_C2 = ',num2str(Snr_C2)]);


