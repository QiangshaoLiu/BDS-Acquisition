clear all;
close all;
%%导频分量的伪码生成，已验证
N = 10243; %legendre序列长度
w = 796;   %B1C导频分量主码PRN=1
p = 7575;  %截取点
n =10230;  %B1C导频分量主码长度
L = zeros(1,N);
W = zeros(1,n);

for k = 1 : N-1 
    for x = 1 : (N-1)/2
    if mod(k,N) == mod(x^2, N)
        L(k) = 1;
    end
    end
end

L = [0 L(1:N-1)];   %生成legendre序列 

%生成B1C导频分量主码
for k = 0 : n-1     
    W(k+1) = mod(sum([L(mod((k+p-1), N)+1), L(mod((k+w+p-1), N)+1)]),2);
end

%生成B1C导频分量子码
N1 = 3607;
w1 = 269;
p1 = 1889;
n1 = 1800;  %B1C导频分量子码长度
L1 = zeros(1,N1);
W1 = zeros(1,n1);

for k = 1 : N1-1 
    for x = 1 : (N1-1)/2
    if mod(k,N1) == mod(x^2, N1)
        L1(k) = 1;
    end
    end
end

L1 = [0 L1(1:N1-1)];   %生成legendre序列


for k = 0 : n1-1     
    W1(k+1) = mod(sum([L1(mod((k+p1-1), N1)+1), L1(mod((k+w1+p1-1), N1)+1)]),2);
end

%%伪码的性能分析
%将Weil码变换成极性码，0表示高电平‘+1’，1表示低电平‘-1’
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
max_cor=max(prn_cor);
%prn_cor_nor = prn_cor/max(prn_cor);
prn_cor_nor = prn_cor/10230;
prn_prn_corr = fftshift(prn_cor_nor);
figure(1)
%plot(prn_prn_corr);title('Weil码自相关');
index=-5115:5114;
plot(index,prn_prn_corr);title('Weil码自相关函数');
%axis([0 10230 -0.05 1.05]);
axis([-5115 5114 -0.05 1.05]);
%生成PRN=2的Weil码用于互相关运算
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
max_cr_cor=max(prn_cor1);
%prn_cor_nor1 = prn_cor1/max(prn_cor1);
prn_prn_corr1 = fftshift(prn_cor1);
figure(2)
%plot(prn_prn_corr1);
plot(index,prn_prn_corr1);title('Weil码互相关函数');
axis([-5115 5114 -500 500]);

figure(3)
subplot(2,1,1);
plot(index,prn_prn_corr);title('Weil码自相关函数');
xlabel('延时/码片');
axis([-5115 5114 -0.05 1.05]);
subplot(2,1,2);
plot(index,prn_prn_corr1);title('Weil码互相关函数');
xlabel('延时/码片');
axis([-5115 5114 -450 450]);


