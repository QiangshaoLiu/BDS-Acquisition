function settings = initSettings()
%% Processing settings ===========================================
settings.msToProcess = 0.1;        %[s]

% Front-End settings
% feID =
% 1 - 毕设自定义
feID = 1;
settings = FrontEndsettings(settings,feID);

%% Constellation characteristics
settings.GNSS_signal = 'BDS_B1C';

% Setting of Integration time [s] for acquisition
switch settings.GNSS_signal
    case ('BDS_B1C')
        settings.T_int = 10e-3; % [s] coherent integration time
        settings.codeFreqBasis      = 1.023e6;      %[Hz]
        % Define number of chips in a code period
        settings.codeLength         = 10230;
        otherwise
        error('err:unkSignal','ERROR: Unknown GNSS signal!');
end

%% Acquisition settings ==========================================
% Band around IF to search for satellite signal. Depends on max Doppler
settings.acqSearchBand      = 10000;       %[Hz]
settings.acqSearchStep      = 250;         %[Hz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 1.2;  %捕获门限值
settings.Non_Coh_Sums       = 5;    %非相干积分
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition
settings.acqSatelliteList   = [1 4 11 15 18 21 24];

%% Constants =====================================================
settings.f0_B1C             = 1575.42e6;      % [Hz] B1C radio-frequency
settings.fm_B1C             = 24.58e6;        % [Hz] B1C medium-frequency
settings.datacode_delay     = 224;            % [chip] 数据通道延迟码相位,必须>1
settings.fd_carr            = 1000;           % [Hz] 载波多普勒频移   
