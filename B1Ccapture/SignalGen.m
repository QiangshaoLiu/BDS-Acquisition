clc
clear all
close all

%%Weil������
N = 10243; %legendre���г���
w = 796;   %B1C��Ƶ��������PRN=1
p = 7575;  %��ȡ��
n =10230;  %B1C��Ƶ�������볤��
L = zeros(1,N);
W = zeros(1,n);

for k = 1 : N-1 
    for x = 1 : (N-1)/2
    if mod(k,N) == mod(x^2, N)
        L(k) = 1;
    end
    end
end

L = [0 L(1:N-1)];   %����legendre���� 

%����B1C��Ƶ��������
for k = 0 : n-1     
    W(k+1) = mod(sum([L(mod((k+p-1), N)+1), L(mod((k+w+p-1), N)+1)]),2);
end

%����B1C��Ƶ��������
N1 = 3607;
w1 = 269;
p1 = 1889;
n1 = 1800;  %B1C��Ƶ�������볤��
L1 = zeros(1,N1);
W1 = zeros(1,n1);

for k = 1 : N1-1 
    for x = 1 : (N1-1)/2
    if mod(k,N1) == mod(x^2, N1)
        L1(k) = 1;
    end
    end
end

L1 = [0 L1(1:N1-1)];   %����legendre����

for k = 0 : n1-1     
    W1(k+1) = mod(sum([L1(mod((k+p1-1), N1)+1), L1(mod((k+w1+p1-1), N1)+1)]),2);
end

%�����
n2 = n * n1;
C = zeros(1, n2);
ii = 1;

for k = 1 : n1
     for i = 1 : n
         C(ii) = mod(sum([W(i),W1(k)]),2); 
         ii = ii + 1;
     end
end

for k = 1 : n2
    if C(k) == 0
        C(k) = 1;
    else
        C(k) = -1;
    end
end

for k = 1 : n
    if W(k) == 0
        W(k) = 1;
    else
        W(k) = -1;
    end
end

%%BOC�ź�����
f_medium = 20*1.023e6;      %��ƵƵ��
f_sample = 70*1.023e6;      %����Ƶ��
f_code = 1.023e6;           %������
f_sc_a = 1.023e6;           %���ݷ������ز�����
f_sc_b = 6.138e6;           %��Ƶ�������ز�����
SNR = -20;                  %�����-20db
sim_t = 1e-1;               %����ʱ��100ms
t = 0 : 1/f_sample : sim_t - 1/f_sample;
t_n = f_code * sim_t;       %����ʱ���ڵ�α������
c_fs = f_sample / f_code;   %ÿ����Ԫ������������

for k = 1 : t_n
    i = mod(k, n2);
    if i == 0
        PseduCode(1+(k-1)*f_sample/f_code : k*f_sample/f_code) = C(n);
    else
        PseduCode(1+(k-1)*f_sample/f_code : k*f_sample/f_code) = C(i);
    end
end

SC_pilot_a = sqrt(29/44)*sign(sin(2*pi*f_sc_a*t));
SC_pilot_b = sqrt(1/11)*sign(sin(2*pi*f_sc_b*t));

S_pilot = PseduCode .*(SC_pilot_a+1i*SC_pilot_b).*(cos(2*pi*f_medium*t)+1i*sin(2*pi*f_medium*t));

figure(1)
plot(t, S_pilot);
title('��Ƶ����ʱ����');

Signal_fft = fft(S_pilot);
Signal_pow = Signal_fft.*conj(Signal_fft);

figure(2)
plot(Signal_pow);
title('��Ƶ����Ƶ�׷���');


