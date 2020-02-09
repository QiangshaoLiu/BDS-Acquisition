function settings = FrontEndsettings(settings, feID)
% INPUT:
% settingsID: front end ID, integer for identifying front-end
% OUTPUT:
% settings: front-end data inside settings structure

settings.snr = -25;       %[db]

switch feID
    case (1) %LSQ
        settings.signalType= 'IF';
        settings.NSample = 1;
        % Intermediate frequency
        settings.IF                 = 24.58e6;  % [Hz]
        % Sampling and code frequencies
        settings.samplingFreq       = 90e6;     % [Hz]
    otherwise
        error('err:NoFrontend','ERROR: No front end!')
end