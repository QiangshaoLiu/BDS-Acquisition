function Subcarr = GenerateSubcarr(flag, settings)

f_sample = settings.samplingFreq;            %采样频率
f_sc_a = settings.codeFreqBasis ;           %数据分量子载波速率
f_sc_b = 6*settings.codeFreqBasis ;           %导频分量子载波速率
sim_t = settings.msToProcess;
t = 0 : 1/f_sample : sim_t - 1/f_sample;

j=sqrt(-1);

switch flag
    case (1)
        Subcarr = sign(sin(2*pi*f_sc_a*t));
    case (2)
        Subcarr = sqrt(29/33)*sign(sin(2*pi*f_sc_a*t))...
            -j*sqrt(4/33)*sign(sin(2*pi*f_sc_b*t));
    otherwise
        error('err:unkBOCsignal','ERROR: Unknown BOC signal!')
end



