%Author:LSQ
%Date:2019/4
%Description: 论文3.2.5章节.

clc;
close all;

set(0,'defaultfigurecolor','w'); %将仿真图背景设置为白色

%利用BOC(alfa,beta)信号功率谱密度的公式绘图
alfa = 15;
beta = 10;
f0 = 1.023e6;
pi = 3.1415926535898;
fs = alfa*f0;
Rc = beta*f0;
n = 2*alfa/beta;
ts = 1/(2*fs);
tc = 1./(beta.*f0);

f = -20./tc:0.0001./tc:20./tc;
num = length(f);

Psd_boc = ts./n.*sinc(f.*ts).^2.*(2*cos(2*pi.*f.*ts)-1).^2; 

x = -4./tc:0.0001./tc:4./tc;
num_x = length(x);
Upper_filter = zeros(1,num_x);
Lower_filter = zeros(1,num_x);

for i = 1:num_x
    if (0.00145./tc<=i)&&(i<=0.0034./tc)
        Upper_filter(i) = -73.48;
    else
        Upper_filter(i) = -110;
    end
end

for i = 1:num_x
    if (0.0044./tc<=i)&&(i<=0.0064./tc)
        Lower_filter(i) = -73.48;
    else
        Lower_filter(i) = -110;
    end
end

figure(1)
plot(f,10*log10(Psd_boc));
set(gca,'xticklabel',{'-40','-30','-20','-10','0','10','20','30','40'});%转化为MHz
hold on;
plot(x,Upper_filter,'r--',x,Lower_filter,'g--');
axis([-4./tc,4./tc,-110,-70]);
xlabel('频率(MHz)');
ylabel('功率谱密度(dBW/Hz)');
legend('BOC(15,10)','上边带滤波','下边带滤波');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
%title('BOC(15,10)信号上下边带滤波');


