%--------------------------------------------------------------------------
%���ڣ�2018��3��
%���ߣ�LSQ
%���ݣ�����B1C�źŵĲ���
%�ź�Ƶ�ʣ�1575.42MHz
%�źŴ���32.736MHz
%��ƵƵ�ʣ�24.58MHz
%����Ƶ�ʣ�90MHz
%�ز�������Ƶ�ƣ�1550Hz
%����ͨ������λ�ӳ٣�224����Ƭ
%��Ƶͨ������λ�ӳ٣�374����Ƭ
%����ȣ�-25dB
%Ƶ������������250Hz
%�������ʱ�䣺10ms
%--------------------------------------------------------------------------
clc;
close all;
disp ('Starting processing...');

settings = initSettings();

% Carrier frequencies of detected signals
acqResults.carrFreq = 0;
% PRN code phases of detected signals
acqResults.codePhase = 0;
% Correlation peak ratios of the detected signals
acqResults.peakMetric = 0;

PRN = 1;

data = GenerateSignal(PRN,settings);
[doppler_est,code_phase,peakMetric]=signal_acquisition_FFT(data,settings,PRN);

acqResults.carrFreq(PRN)  = doppler_est;
acqResults.codePhase(PRN) = code_phase;
acqResults.peakMetric(PRN) = peakMetric;
acqResults.acqThreshold = settings.acqThreshold;

disp ('Processing done!');

