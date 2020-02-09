%-------------------------------------------------------------------------
% ������Ϊ����
% ���ߣ�LSQ
% ���ڣ�2016��10��29��
%----------------------����˵��-------------------------------------------
% ������ʵ��GPS L1�źŵ������������񲢸���������,��GXQ�ĳ���Ϊ������д��
% �޸Ķ�Ӧ�������޸Ķ�����Ƶ�ƣ���Ƶ�����ʣ���ʼ����λƫ�Ƶ�ʺ�����λ��������
%----------------------�������-------------------------------------------
% L1��Ƶ1575.42M��ģ����Ƶ��1.5M��C/A�����ʣ�1.023M��������Ƶ�ƣ�1K��
% �����ռ��ٶȣ�0Hz/s���������볤�ȣ�2ms����ʼ�ز���λƫ�0rad��
% ��ʼ����λ����481��Ƭ����Ƶ����Ƶ�ʣ�5MHz�������źŹ��ʣ�-130dBm��
% �����ź�����ȣ�-20dB��Ƶ��������Χ��1.5M-2K~1.5M+2K��Ƶ������������100Hz��
% ����λ����������0.5��Ƭ����ɻ���ʱ�䣺1ms������ɻ�����Ŀ��1
%-------------------------------------------------------------------------
fm=1.5e6;                                       %��Ƶ
fs=5e6;                                         %����Ƶ��
fd=1e3;                                         %������Ƶ��
fbin=100;                                       %Ƶ����������
sear_n=2;                                       %Ƶ����������Ƶ��������յ�Ƶ�Ƶı�����ϵ����Ӧ������ΧΪf0-sear_n*fd~f0+sear_n*fd
delay=480;                                      %�ӳٵ���Ƭ��������Ϊ�ӳ�480��Ƭ�������յ���CA��Ƭ�ӵ�481��Ƭ��ʼ
SNR=-20;                                        %��Ƶ����������-20db
snr=10^(SNR/10);                                %����ȣ���ͨ�ı�����ϵ
sim_t=2;                                        %����ʱ��2ms
t=0:1e3/fs:sim_t-1e3/fs;                        %�����������ܳ��ȣ�2ms����10000��������
t_n=length(t)/fs;                               %���������ݳ��ȣ���λ���룩
L=length(t);                                    %�ܲ�������
l=fs*1e-3;                                      %1ms��Ӧ�ĵ������ڴ�Ϊ5000
PRN=5;                                          %���Ǳ��
Pfe=1e-6;                                       %�龯��
noi=randn(1,L);                                 %������˹������
Power_noi=noi*noi'/L;                           %��˹���������� noi'Ϊnoi��ת��
Vt=sqrt(Power_noi)*sqrt(-2*log(Pfe));           %�������

%����C/A�뷢�����������������ݳ���һ�µ�α������
[CAcode_delay]=CAgen(fs,delay,sim_t);       %�������Ǳ��Ϊ5��α������CAcode�Լ���fs��������2ms��CA�����ֵ

%��������λ
Data_pre=ones(1,l/2);                       %1ms��Ӧ�ĵ������Ĳ�����
Data=[Data_pre,repmat(-1*Data_pre,1,3)];    

%�������ǵ����ź�
Power_xifs=Power_noi*snr;
A=sqrt(2*Power_xifs);
xif_s=A*Data.*CAcode_delay.*sin(2*pi*(fm+fd)*t);        %��Ƶ���ǵ����ź�
figure(1);subplot(321);plot(t,xif_s);
grid on;title('��Ƶ���ǵ����ź�ʱ��ͼ');xlabel('t');ylabel('xif_s(t)');
N=2^(nextpow2(L));
M=(0:N-1)*fs/N;
X=fft(xif_s,N);
figure(1);subplot(322);plot(M,abs(X));
grid on;title('��Ƶ���ǵ����ź�Ƶ��ͼ');xlabel('f');ylabel('xif_s(f)');

%��Ƶ�����ź�
xif=xif_s+noi;                                  %��Ƶ�ź�����Ϊ�����źź�����
figure(1);subplot(323);plot(t,noi);
grid on;title('��˹������ʱ��ͼ');xlabel('t');ylabel('noise');
figure(1);subplot(325);plot(t,xif);
grid on;title('��Ƶ�����ź�ʱ��ͼ');xlabel('t');ylabel('xif(t)');
X_noise=fft(xif,N);noise=fft(noi,N);
figure(1);subplot(324);plot(M,abs(noise));
grid on;title('��˹������Ƶ��ͼ');xlabel('f');ylabel('noise');
figure(1);subplot(326);plot(M,abs(X_noise));
grid on;title('��Ƶ�����ź�Ƶ��ͼ');xlabel('f');ylabel('xif(f)');

%�������ջ����ػָ���Ƶ�ź�
loop_t=0.5;                                     %��������λ��������λ����Ƭ
loop_chip=0:loop_t:1023-loop_t;                 %����λһ����������Ƭ                      
loop_n=length(loop_chip);                       %����λ�����ܵ�Ԫ������0��Ƭ��ʼ����1023��Ƭ����
f0=fm;                                          %���ػָ��ز������ĸ�����Ƶ
t1=0:1e3/fs:1-1e3/fs;
t2=1:1e3/fs:2-1e3/fs;

%ѭ��������������λ��Ƶ��
CA_loop_1=zeros(1,l);                           %���ѭ��ƽ�ƺ��α��
CAcatch=zeros(sear_n*2*fd/fbin+1,loop_n);       %���ڴ������������Ԫ��󲶻�������ֵ
for m=-sear_n*fd/fbin:sear_n*fd/fbin            %m���Ʊ������е�����Ƶ�ʵ�
    for n=loop_chip                             %n���Ʊ������е���������λ��
        [CAcode_delay]=CAgen(fs,n,sim_t);       %����CA���������һ�β������ӳ�delay��Ƭ����������Ƭ���ӳ��ɲ���n����
        CA_loop1=CAcode_delay(1:l);
        CA_loop2=CAcode_delay(l+1:L);
        xli1=CA_loop1.*sin(2*pi*(f0+m*fbin)*t1);
        xli2=CA_loop2.*sin(2*pi*(f0+m*fbin)*t2);
        xlq1=CA_loop1.*cos(2*pi*(f0+m*fbin)*t1);
        xlq2=CA_loop2.*cos(2*pi*(f0+m*fbin)*t2);
        xcorri1=sum(xif(1:l).*xli1)/l;                   %��I֧·���ֵ�����е�1ms��������
        xcorrq1=sum(xif(1:l).*xlq1)/l;                   %��Q֧·���ֵ�����е�1ms��������
        xcorri2=sum(xif(l+1:L).*xli2)/l;                 %��I֧·���ֵ�����е�2ms��������
        xcorrq2=sum(xif(l+1:L).*xlq2)/l;                 %��Q֧·���ֵ�����е�2ms��������
        xcorri=max(xcorri1,xcorri2);         %���������λ����ֱ����������Ļ������㣬��������������1����������ǲ���������ģ�������ֵ��������λӰ��
        xcorrq=max(xcorrq1,xcorrq2);         %ͬ��
        xncor=abs(xcorri+1i*xcorrq);                   %���з���ɻ�����Ϊ1�ķ���ɻ��֣�������ģ����
        CAcatch(m+sear_n*fd/fbin+1,n/loop_t+1)=xncor;  %�����ճ���Ϊ1�ķ���ػ���ֵ������Ӧ��������Ԫ��Ϊ����ֵ�ļ������     
    end
end
figure(2);plot(t1,xli1);
grid on;title('I֧·�ź�ʱ��ͼ');xlabel('t');ylabel('xlo(t)');
figure(3);surf(CAcatch);title('��ά����ͼ');
xlabel('��������λƫ��');ylabel('����Ƶ��ƫ��');zlabel('����ֵ');

%Ѱ�Ҳ������в���Ԫ�����ֵ�±����������������Ƶ�ƺ�����λƫ��
[cat_fidex,cat_caidex]=find(CAcatch==max(max(CAcatch)));
cat_f=(cat_fidex-(sear_n*fd/fbin+1))*fbin;          %�������Ƶ��
cat_ca=(cat_caidex-1)*loop_t;                       %������λƫ��

%��������ʾ
disp(['Catch Doppler frequency shift is:fd=',num2str(cat_f),'Hz']);
disp(['Cacth CA Code shift is:CAcode delay is=',num2str(cat_ca),'chip']);

