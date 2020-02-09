function [Code, Rc] = GenerateRecvCode(PRN, GNSS_signal, datacode_delay)
%%clc;
%close all;
%GNSS_signal = 'BDS_B1C';
%PRN = 1;
%datacode_delay = 4;

switch GNSS_signal
    case ('BDS_B1C')
        Code = Legen(PRN);         % BDS local code generator 
        Rc = 1.023e6;              % BDS data code rate.
        
        %模拟卫星信号，设定码相位时延
        Code = [Code(datacode_delay :10230) Code(1:datacode_delay-1)];
    otherwise
        error('err:unkSignal','ERROR: Unknown GNSS signal!')
end




