%���γ���������֤�����ģ���е�FFT IP���Ƿ��������
%data = [1,2,3,4,5,6,7,8];
%fid = fopen('D:\My_project\project_5\data.txt','wt');
%fprintf(fid,'%x\n',data);
%fclose(fid);

f_sample = 8.192e6;             %����Ƶ��
T_process = 1e-3;                 %����ʱ��
f_sc_a = 1.023e6;           %���ݷ������ز�����

t = 0 : 1/f_sample : T_process - 1/f_sample;

%Subcarr1 = zeros(1,8192);
%for k=1:8192
Subcarr1 = sin(2*pi*f_sc_a*t);
%if mod(k,2)==0;
%    Subcarr1(k) = 1;
%end
%end
Rc = 1.023e6;                      %����������
prn_p = generatecode(2);
index_code = mod(floor(Rc*t),10230)+1;
prn_local = prn_p(index_code);

q = quantizer('single');   %ת��ΪIEEE-754�����ȸ�������������ʽ
%temp = num2bin(q, Subcarr1);
%data = bin2dec(temp);      %ת��Ϊʮ����������ʮ�����Ƹ�ʽд��txt�ļ�
value = abs(fft(prn_local));
value1 = abs(fft(prn_local).*fft(prn_local));
tempI = num2bin(q, value);
tempQ = num2bin(q, value1);
dataI = bin2dec(tempI);
dataQ = bin2dec(tempQ);

%fid = fopen('D:\My_project\project_5\data.txt','wt');
%fprintf(fid,'%x\n',data);
%fclose(fid);
  fidI = fopen('D:\My_project\project_5\codeI.txt','wt');
  fprintf(fidI,'%x\n',dataI);
  fclose(fidI);
  fidQ = fopen('D:\My_project\project_5\codeQ.txt','wt');
  fprintf(fidQ,'%x\n',dataQ);
  fclose(fidQ);

figure(1)
plot(Subcarr1);title('�ز�ʱ����');

fft_code = fft(prn_local);
fft_mul = fft_code.*fft_code;

figure(2)
plot(fft_code);



