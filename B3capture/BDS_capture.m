%-------------------------------------------------------------------------
% 化腐朽为神奇
% 作者：LSQ
% 日期：2016年10月29日
%----------------------程序说明-------------------------------------------
% 本程序实现GPS L1信号的线性搜索捕获并给出捕获结果,以GXQ的程序为基础编写。
% 修改对应参数可修改多普勒频移，中频采样率，初始码相位偏差，频率和码相位搜索步长
%----------------------程序参数-------------------------------------------
% L1载频1575.42M；模拟中频：1.5M；C/A码速率：1.023M；多普勒频移：1K；
% 多普勒加速度：0Hz/s；数据输入长度：2ms；初始载波相位偏差：0rad；
% 初始码相位：第481码片；中频采样频率：5MHz；接收信号功率：-130dBm；
% 接收信号信噪比：-20dB；频率搜索范围：1.5M-2K~1.5M+2K；频率搜索步长：100Hz；
% 码相位搜索步长：0.5码片；相干积分时间：1ms；非相干积分数目：1
%-------------------------------------------------------------------------
fm=1.5e6;                                       %中频
fs=5e6;                                         %采样频率
fd=1e3;                                         %多普勒频移
fbin=100;                                       %频域搜索步长
sear_n=2;                                       %频率搜索单边频宽与多普勒的频移的倍数关系，对应搜索范围为f0-sear_n*fd~f0+sear_n*fd
delay=480;                                      %延迟的码片数，设置为延迟480码片，即接收到的CA码片从第481码片开始
SNR=-20;                                        %中频输入的信噪比-20db
snr=10^(SNR/10);                                %信噪比，普通的倍数关系
sim_t=2;                                        %仿真时间2ms
t=0:1e3/fs:sim_t-1e3/fs;                        %产生的数据总长度，2ms共计10000个采样点
t_n=length(t)/fs;                               %产生的数据长度（单位毫秒）
L=length(t);                                    %总采样点数
l=fs*1e-3;                                      %1ms对应的点数，在此为5000
PRN=5;                                          %卫星编号
Pfe=1e-6;                                       %虚警率
noi=randn(1,L);                                 %产生高斯白噪声
Power_noi=noi*noi'/L;                           %高斯白噪声功率 noi'为noi的转置
Vt=sqrt(Power_noi)*sqrt(-2*log(Pfe));           %检测门限

%调用C/A码发生器产生与输入数据长度一致的伪码数据
[CAcode_delay]=CAgen(fs,delay,sim_t);       %产生卫星编号为5的伪卫星码CAcode以及在fs采样率下2ms的CA码采样值

%产生数据位
Data_pre=ones(1,l/2);                       %1ms对应的导航电文采样点
Data=[Data_pre,repmat(-1*Data_pre,1,3)];    

%产生卫星调制信号
Power_xifs=Power_noi*snr;
A=sqrt(2*Power_xifs);
xif_s=A*Data.*CAcode_delay.*sin(2*pi*(fm+fd)*t);        %中频卫星调制信号
figure(1);subplot(321);plot(t,xif_s);
grid on;title('中频卫星调制信号时域图');xlabel('t');ylabel('xif_s(t)');
N=2^(nextpow2(L));
M=(0:N-1)*fs/N;
X=fft(xif_s,N);
figure(1);subplot(322);plot(M,abs(X));
grid on;title('中频卫星调制信号频域图');xlabel('f');ylabel('xif_s(f)');

%中频接收信号
xif=xif_s+noi;                                  %中频信号输入为卫星信号和噪声
figure(1);subplot(323);plot(t,noi);
grid on;title('高斯白噪声时域图');xlabel('t');ylabel('noise');
figure(1);subplot(325);plot(t,xif);
grid on;title('中频卫星信号时域图');xlabel('t');ylabel('xif(t)');
X_noise=fft(xif,N);noise=fft(noi,N);
figure(1);subplot(324);plot(M,abs(noise));
grid on;title('高斯白噪声频域图');xlabel('f');ylabel('noise');
figure(1);subplot(326);plot(M,abs(X_noise));
grid on;title('中频卫星信号频域图');xlabel('f');ylabel('xif(f)');

%产生接收机本地恢复载频信号
loop_t=0.5;                                     %搜索码相位步长，单位：码片
loop_chip=0:loop_t:1023-loop_t;                 %码相位一次搜索的码片                      
loop_n=length(loop_chip);                       %码相位搜索总单元数，从0码片开始搜起，1023码片结束
f0=fm;                                          %本地恢复载波产生的复制中频
t1=0:1e3/fs:1-1e3/fs;
t2=1:1e3/fs:2-1e3/fs;

%循环遍历所有码相位和频率
CA_loop_1=zeros(1,l);                           %存放循环平移后的伪码
CAcatch=zeros(sear_n*2*fd/fbin+1,loop_n);       %用于存放所有搜索单元最后捕获运算后的值
for m=-sear_n*fd/fbin:sear_n*fd/fbin            %m控制遍历所有的搜索频率点
    for n=loop_chip                             %n控制遍历所有的搜索码相位点
        [CAcode_delay]=CAgen(fs,n,sim_t);       %调用CA码产生函数一次产生从延迟delay码片整数倍的码片，延迟由参数n控制
        CA_loop1=CAcode_delay(1:l);
        CA_loop2=CAcode_delay(l+1:L);
        xli1=CA_loop1.*sin(2*pi*(f0+m*fbin)*t1);
        xli2=CA_loop2.*sin(2*pi*(f0+m*fbin)*t2);
        xlq1=CA_loop1.*cos(2*pi*(f0+m*fbin)*t1);
        xlq2=CA_loop2.*cos(2*pi*(f0+m*fbin)*t2);
        xcorri1=sum(xif(1:l).*xli1)/l;                   %求I支路相关值并进行第1ms积分运算
        xcorrq1=sum(xif(1:l).*xlq1)/l;                   %求Q支路相关值并进行第1ms积分运算
        xcorri2=sum(xif(l+1:L).*xli2)/l;                 %求I支路相关值并进行第2ms积分运算
        xcorrq2=sum(xif(l+1:L).*xlq2)/l;                 %求Q支路相关值并进行第2ms积分运算
        xcorri=max(xcorri1,xcorri2);         %解决跨数据位问题分别进行两毫秒的积分运算，两毫秒中至少有1毫秒的数据是不发生跳变的，即积分值不受数据位影响
        xcorrq=max(xcorrq1,xcorrq2);         %同上
        xncor=abs(xcorri+1i*xcorrq);                   %进行非相干积分数为1的非相干积分，复数求模运算
        CAcatch(m+sear_n*fd/fbin+1,n/loop_t+1)=xncor;  %将最终长度为1的非相关积分值赋予相应的搜索单元作为捕获值的检测输入     
    end
end
figure(2);plot(t1,xli1);
grid on;title('I支路信号时域图');xlabel('t');ylabel('xlo(t)');
figure(3);surf(CAcatch);title('三维曲面图');
xlabel('捕获码相位偏差');ylabel('捕获频率偏差');zlabel('捕获值');

%寻找捕获结果中捕获单元的最大值下标索引并计算多普勒频移和码相位偏差
[cat_fidex,cat_caidex]=find(CAcatch==max(max(CAcatch)));
cat_f=(cat_fidex-(sear_n*fd/fbin+1))*fbin;          %算多普勒频移
cat_ca=(cat_caidex-1)*loop_t;                       %算码相位偏差

%捕获结果显示
disp(['Catch Doppler frequency shift is:fd=',num2str(cat_f),'Hz']);
disp(['Cacth CA Code shift is:CAcode delay is=',num2str(cat_ca),'chip']);

