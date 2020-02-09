%本程序利用BOC(alfa,beta)信号功率谱密度的公式绘图

clc;
close all;

alfa = 15;
beta = 10;
f0 = 1.023e6;
pi = 3.1415926535898;
fs = alfa*f0;
Rc = beta*f0;
n = 2*alfa/beta;
ts = 1/(2*fs);

f = -50*f0:0.0001*f0:50*f0;
num = length(f);

Psd_boc = 1/n/ts*(sin(pi*ts.*f).*cos(pi*ts*n.*f)./(pi.*f.*cos(pi*ts.*f))).^2;
Amp_max = max(Psd_boc(1:num/2)); 

x = -40 : 1 :40;
num_x = length(x);

Upper_filter = zeros(1,num_x);
Lower_filter = zeros(1,num_x);

for i = 1:num_x
    if (47<=i)&&(i<=65)
        Upper_filter(i) = 1;
    end
end

for i = 1:num_x
    if (17<=i)&&(i<=35)
        Lower_filter(i) = 1;
    end
end

figure(1)
plot(f./10^6,Psd_boc/Amp_max);
hold on;
plot(x,Upper_filter,'r--',x,Lower_filter,'g--');
axis([-40,40,0,1.1]);
xlabel('频率  MHz');
ylabel('归一化功率谱密度 dBW/Hz');
legend('BOC(15,10)','上边带滤波','下边带滤波');
title('BOC(15,10)信号上下边带滤波');
