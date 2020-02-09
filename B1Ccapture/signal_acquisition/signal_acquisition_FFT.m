% FFT Acquisition in Time domain
function [doppler,code_phase,peakMetric]=signal_acquisition_FFT(data,settings,prn)
T_int=settings.T_int;
sampling_freq=settings.samplingFreq;
IF=settings.IF;

j=sqrt(-1);

Non_Coh_Sums = settings.Non_Coh_Sums;

Dopplerstep  = settings.acqSearchStep;
DopplerRange = settings.acqSearchBand;

FD_vect= -DopplerRange:Dopplerstep:DopplerRange;
Doppler_CENTER = 0;
FD_vect = FD_vect + Doppler_CENTER;

acq_metric = settings.acqThreshold;   % Acquisition metric

% Local code generation
[Loc, Rc] = GenerateLocCode(prn,settings.GNSS_signal);

% Local BOC generation
BOC = GenerateSubcarr(1, settings);

%f_sample = settings.samplingFreq;            %采样频率
%sim_t = settings.msToProcess;
%t = 0 : 1/f_sample : sim_t - 1/f_sample;
%fm = settings.fm_B1C;

%Loc = Loc(mod(floor(t*Rc),10230)+1);

%for f_search = fm-DopplerRange:Dopplerstep:fm+DopplerRange
 %   signal1 = fft(data.*sin(2*pi*f_search*t));   %经过载波剥离的信号
   % signal2 = Loc.*GenerateSubcarr(1,settings);  %本地BOC(1,1)信号
 %   value = abs(ifft(signal1.*conj(fft(Loc)))); 
%end

CodeLen = length(Loc);
num_samples = sampling_freq/Rc;         % Number of samples per chip.
Fs = sampling_freq;
fif = rem(IF,Fs);
data = data(1:1:end);

N = floor(Fs*T_int);   % Number of samples for each coherent integration
N_tot = N * (Non_Coh_Sums); % Total number of processed samples

% FFT acquisition
idx = 1;

n = 0:N_tot-1;
ind_cod = mod(floor(n*Rc/Fs),CodeLen)+1;
SigLOC_tot = Loc(ind_cod);

C = zeros(length(FD_vect),N);

for ind_FD= 1:length(FD_vect)
    FD = FD_vect(ind_FD);
    corr = zeros(1,N)+j*zeros(1,N);
    
    k = 0:N-1;
    argx = 2*pi*(fif+FD)/Fs;
    carrI = cos(argx*k);
    carrQ = sin(argx*k);
    
    for M = 0:(Non_Coh_Sums-1)
        SigLOC = SigLOC_tot(N*M+1:N*M+N);
        SigLOCFFT = conj(fft(SigLOC,N));      
        
        BocLOCFFT = conj(fft(BOC,N));
        
        SigIN = data(N*M+1:N*M+N);
        I = SigIN.*carrI;
        Q = SigIN.*carrQ;
        
        SigINIQ = I+j*Q;
        
        corr = corr + abs(ifft(fft(SigINIQ,N).*(BocLOCFFT))) - ...
            abs(ifft(fft(SigINIQ,N).*(SigLOCFFT)));
        
        figure(2)
        
    end
    
    C(idx,:) = corr;
    idx = idx+1;
end

%Find the main peak in the correlation floor and the corresponding frequency bin index
[bb, ind_mixf] = max(max(C'));
[bb, ind_mixc] = max(max(C));

if (ind_mixc < ceil(num_samples)),
    vect_search_peak = [zeros(1,2*ceil(num_samples)), ...
        C(ind_mixf,(2*ceil(num_samples)):end)];
elseif (ind_mixc > (length(C(ind_mixf,:) - 2*ceil(num_samples))))
    vect_search_peak = [C(ind_mixf,1:(end-2*ceil(num_samples)):end),...
        zeros(1,2*ceil(num_samples))];
else
    vect_search_peak = [C(ind_mixf,1:(ind_mixc-ceil(num_samples))),...
        zeros(1,2*ceil(num_samples)-1),...
        C(ind_mixf,(ind_mixc+ceil(num_samples)):end)];
end

%Find the second highest peak in the correlation floor
second_peak = max(vect_search_peak);

%compare the acquisition metric to a predefined threshold
peakMetric=bb/second_peak;

code_phase = floor((ind_mixc-1)*10230/N);
doppler =(ind_mixf-1)*Dopplerstep - 10000;   %[HZ]



