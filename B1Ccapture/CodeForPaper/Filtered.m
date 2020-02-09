%Author:LSQ
%Date:2019/4
%Description: ����3.2.6�½�.

clc;
close all;

set(0,'defaultfigurecolor','w'); %������ͼ��������Ϊ��ɫ

%��֤Filtered�����㷨����ȷ��
legendrelength = 10243;
WeilCodelength = 10230;
legendre = zeros(1, legendrelength);
WeilCode = zeros(1, WeilCodelength);
WeilCode_p = zeros(1, WeilCodelength);
phase_diff = 2678;
intercept_pot = 699;
phase_diff_p = 796;   %PRN=1�ĵ�Ƶ����������λ��
intercept_pot_p = 7575;  %PRN=1�ĵ�Ƶ���������ȡ��

for k = 1 : legendrelength-1 
    for x = 1 : (legendrelength-1)/2
        if mod(k,legendrelength) == mod(x^2, legendrelength)
            legendre(k) = 1;
        end
    end
end

legendre = [0 legendre(1:legendrelength-1)];   %����legendre���� 

%����B1C���ݷ�������
for k = 0 : WeilCodelength-1     
    WeilCode(k+1) = mod(sum([legendre(mod((k+intercept_pot-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff+intercept_pot-1), legendrelength)+1)]),2);
end

%��Weil��任�ɼ����룬0��ʾ�ߵ�ƽ��+1����1��ʾ�͵�ƽ��-1��
for k = 1 : WeilCodelength
    if WeilCode(k) == 0
       WeilCode(k) = 1;
    else
       WeilCode(k) = -1;
    end
end

%����B1C��Ƶ��������
for k = 0 : WeilCodelength-1     
    WeilCode_p(k+1) = mod(sum([legendre(mod((k+intercept_pot_p-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff_p+intercept_pot_p-1), legendrelength)+1)]),2);
end

%��Weil��任�ɼ����룬0��ʾ�ߵ�ƽ��+1����1��ʾ�͵�ƽ��-1��
for k = 1 : WeilCodelength
    if WeilCode_p(k) == 0
       WeilCode_p(k) = 1;
    else
       WeilCode_p(k) = -1;
    end
end

f_sample = 40*1.023e6;             %����Ƶ��
f_sc_a = 1.023e6 ;                 %���ݷ������ز�����
f_sc_b = 6*1.023e6 ;               %��Ƶ�������ز�����
T_process = 10e-3;                 %����ʱ��

t = 0 : 1/f_sample : T_process - 1/f_sample;

j=sqrt(-1);

Subcarr1 = sign(sin(2*pi*f_sc_a*t));
Subcarr2 = sign(sin(2*pi*f_sc_b*t));

Rc = 1.023e6;       %����������

% FFT acquisition
idx = 1;
n=0:length(t)-1;
ind_cod = mod(floor(n*Rc/f_sample),WeilCodelength)+1;
SigLOC_tot_d = WeilCode(ind_cod);
SigLOC_tot_p = WeilCode_p(ind_cod);

num_sample = floor(f_sample/Rc);

[BOCPRN_corr lag] = xcorr(Subcarr1.*SigLOC_tot_d, SigLOC_tot_d, 'coeff');  %BOC�ź���PRN��Ļ���غ���
BOC_corr = xcorr(Subcarr1.*SigLOC_tot_d, 'coeff');    %BOC������غ���
%������߶ȱ任
index = lag / num_sample;

%����filter�ķ������õ����ع���غ���
%���ݷ���
PRN_num =   floor(1000/2000*1/Rc*f_sample);                               %��ǰ/�ͺ�����Ƭ����Ӧ�Ĳ�������
zero_filling = ones(1, PRN_num-1);
a_num = floor(1150/2000*1/Rc*f_sample); 
a_filling = ones(1, a_num-1);
PRN_lag = [SigLOC_tot_d(PRN_num:end) zero_filling];        %�ͺ�����Ƭ��PRN
PRN_advance = [a_filling SigLOC_tot_d(1:end-a_num+1)];%��ǰ�����Ƭ��PRN
BOCPRN_corr_filter_a = xcorr(Subcarr1.*SigLOC_tot_d, PRN_advance, 'coeff');  %BOC�ź��볬ǰPRN��Ļ���غ���
BOCPRN_corr_filter_l = xcorr(Subcarr1.*SigLOC_tot_d, PRN_lag, 'coeff');  %BOC�ź��볬ǰPRN��Ļ���غ���
filter = abs(BOC_corr-0.5*(BOCPRN_corr_filter_a-BOCPRN_corr_filter_l)).^2;
b_max = max(filter);

%��Ƶ����
PRN_p_lag = [SigLOC_tot_p(PRN_num:end) zero_filling];        %�ͺ�����Ƭ��PRN
PRN_p_advance = [a_filling SigLOC_tot_p(1:end-a_num+1)];    %��ǰ�����Ƭ��PRN
BOCPRN_corr_p_filter_a = xcorr(sqrt(1/11)*Subcarr2.*SigLOC_tot_p+j*sqrt(29/44)*Subcarr1.*SigLOC_tot_p, PRN_p_advance, 'coeff');  %BOC�ź���PRN��Ļ���غ���
BOCPRN_corr_p_filter_l = xcorr(sqrt(1/11)*Subcarr2.*SigLOC_tot_p+j*sqrt(29/44)*Subcarr1.*SigLOC_tot_p, PRN_p_lag, 'coeff');  %BOC�ź��볬ǰPRN��Ļ���غ���
filter_p = abs(BOCPRN_corr_p_filter_a)+abs(BOCPRN_corr_p_filter_l) - abs(BOCPRN_corr_p_filter_a+BOCPRN_corr_p_filter_l);
r_max = max(filter_p);

figure(1)
plot(index,sqrt(BOC_corr.^2), 'b',index, filter, 'r');
xlabel('��λ��ʱ(��Ƭ)');
ylabel('���ֵ');
axis([-1.5 1.5 -0.5 2.3]);
legend('BOC(1,1)������غ���','Filtered��ز����㷨');
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

