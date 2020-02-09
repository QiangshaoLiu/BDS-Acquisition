function [Code, Rc] = GenerateLocCode(PRN,GNSS_signal)
%
% GenerateLocCode: this function generates the local code for the
% acquisition, depending on the chosen GNSS signal (band & channel) and PRN. 
% It can be customized in order to call user refined code genrators and
% return the code rate (Rc).

%global GNSS_signal;
%clc;
%close all;
%GNSS_signal = 'BDS_B1C';
switch GNSS_signal
    case ('BDS_B1C')
        %Code = Legen(1);         % BDS local code generator
        Code = Legen(PRN);         % BDS local code generator 
        Rc = 1.023e6;              % BDS data code rate.
    otherwise
        error('err:unkSignal','ERROR: Unknown GNSS signal!')
end
