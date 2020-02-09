function signal = GenerateSignal(PRN,settings)

% Recv code generation
[RecvCode, RecvRc] = GenerateRecvCode(PRN, settings.GNSS_signal,settings.datacode_delay);

Signal_flag = 1;   %1:数据分量  2:导频分量
Subcarr_d = GenerateSubcarr(Signal_flag,settings);

f_sample = settings.samplingFreq;            %采样频率
sim_t = settings.msToProcess;
f_code = RecvRc;    %码率
t = 0 : 1/f_sample : sim_t - 1/f_sample;
fm = settings.fm_B1C;
fd = settings.fd_carr;

CodeSample = RecvCode(mod(floor(t*f_code),10230)+1);%时间t内经过采样的weil码
Subcarr_data = CodeSample.*Subcarr_d;            %BOC(1,1)

signal = sin(2*pi*(fm+fd)*t).*Subcarr_data;

%SNR = settings.snr;

%signal = awgn(signal,SNR,'measured');
figure(1)
signal_fft = abs(fft(signal));
plot(signal_fft);
