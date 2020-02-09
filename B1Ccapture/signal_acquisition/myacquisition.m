%--------------------------------------------------------------------------
%日期：2018年3月
%作者：LSQ
%内容：北斗B1C信号的捕获
%信号频率：1575.42MHz
%信号带宽：32.736MHz
%中频频率：24.58MHz
%采样频率：90MHz
%载波多普勒频移：1550Hz
%数据通道码相位延迟：224个码片
%导频通道码相位延迟：374个码片
%信噪比：-25dB
%频率搜索步进：250Hz
%捕获积分时间：10ms
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

