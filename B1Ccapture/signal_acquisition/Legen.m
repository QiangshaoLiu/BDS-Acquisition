function WeilCode = Legen(svnum)
%clc;
%close all;
legendrelength = 10243;
WeilCodelength = 10230;
legendre = zeros(1, legendrelength);
WeilCode = zeros(1, WeilCodelength);
%svnum =1;
switch svnum
    case (1)
        phase_diff = 2678;
        intercept_pot = 699;
    otherwise
        error('err:noCode','ERROR: Unknown legendre code!')
end

for k = 1 : legendrelength-1 
    for x = 1 : (legendrelength-1)/2
        if mod(k,legendrelength) == mod(x^2, legendrelength)
            legendre(k) = 1;
        end
    end
end

legendre = [0 legendre(1:legendrelength-1)];   %生成legendre序列 

%生成B1C数据分量主码
for k = 0 : WeilCodelength-1     
    WeilCode(k+1) = mod(sum([legendre(mod((k+intercept_pot-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff+intercept_pot-1), legendrelength)+1)]),2);
end

%将Weil码变换成极性码，0表示高电平‘+1’，1表示低电平‘-1’
for k = 1 : WeilCodelength
    if WeilCode(k) == 0
       WeilCode(k) = 1;
    else
       WeilCode(k) = -1;
    end
end

