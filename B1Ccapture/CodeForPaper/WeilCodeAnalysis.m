%Author: LSQ
%Date: 2019/4
%Description: ����2.1.1�½�.

clc;
close all;

set(0,'defaultfigurecolor','w'); %������ͼ��������Ϊ��ɫ

%%��Ƶ������α������
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

%%α������ܷ���
%��Weil��任�ɼ����룬0��ʾ�ߵ�ƽ��+1����1��ʾ�͵�ƽ��-1��
for k = 1 : n
    if W(k) == 0
       W(k) = 1;
    else
       W(k) = -1;
    end
end

prn_fft = fft(W);
prn_pow = prn_fft.*conj(prn_fft);
prn_cor = ifft(prn_pow);
prn_prn_corr = fftshift(prn_cor);

index=-5115:5114;

%����PRN=2��Weil�����ڻ��������
L2 = zeros(1,N);
W2 = zeros(1,n);
for k = 1 : N-1 
    for x = 1 : (N-1)/2
    if mod(k,N) == mod(x^2, N)
        L2(k) = 1;
    end
    end
end

L2 = [0 L2(1:N-1)];    

for k = 0 : n-1     
    W2(k+1) = mod(sum([L2(mod((k+2369-1), N)+1), L2(mod((k+156+2369-1), N)+1)]),2);
end

for k = 1 : n
    if W2(k) == 0
       W2(k) = 1;
    else
       W2(k) = -1;
    end
end

prn_fft1 = fft(W2);
prn_pow1 = prn_fft.*conj(prn_fft1);
prn_cor1 = ifft(prn_pow1);
prn_prn_corr1 = fftshift(prn_cor1);

figure(1)
subplot(2,1,1);
plot(index,prn_prn_corr);title('Weil������غ���');
xlabel('��λ��ʱ(��Ƭ)');
ylabel('���ֵ');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([-5115 5114 -1000 10300]);
subplot(2,1,2);
plot(index,prn_prn_corr1);title('Weil�뻥��غ���');
xlabel('��λ��ʱ(��Ƭ)');
ylabel('���ֵ');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([-5115 5114 -450 450]);

figure(2)
plot(index,prn_prn_corr);
xlabel('��λ��ʱ(��Ƭ)');
ylabel('���ֵ');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([-5115 5114 -1000 10500]);

figure(3)
plot(index,prn_prn_corr1);
xlabel('��λ��ʱ(��Ƭ)');
ylabel('���ֵ');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([-5115 5114 -450 450]);


