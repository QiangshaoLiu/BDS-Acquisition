function [ WeilCode_p ] = generatecode( PRN )
%GENERATECODE Summary of this function goes here
%   Detailed explanation goes here
%产生PRN号为1至10的主码序列

switch PRN
    case (1)
        phase_diff_p = 796;
        intercept_pot_p = 7575;
    case (2)
        phase_diff_p = 156;
        intercept_pot_p = 2369;
    case (3)
        phase_diff_p = 4198;
        intercept_pot_p = 5688;
    case (4)
        phase_diff_p = 3941;
        intercept_pot_p = 539;
    case (5)
        phase_diff_p = 1374;
        intercept_pot_p = 2270;
    case (6)
        phase_diff_p = 1338;
        intercept_pot_p = 7306;
    case (7)
        phase_diff_p = 1833;
        intercept_pot_p = 6457;
    case (8)
        phase_diff_p = 2521;
        intercept_pot_p = 6254;
    case(9)
        phase_diff_p = 3175;
        intercept_pot_p = 5644;
    case (10)
        phase_diff_p = 168;
        intercept_pot_p = 7119;
    otherwise
        phase_diff_p = 796;
        intercept_pot_p = 7575;    %超出范围，默认PRN=1
end 

legendrelength = 10243;
WeilCodelength = 10230;
legendre = zeros(1, legendrelength);
WeilCode_p = zeros(1, WeilCodelength);

for k = 1 : legendrelength-1 
    for x = 1 : (legendrelength-1)/2
        if mod(k,legendrelength) == mod(x^2, legendrelength)
            legendre(k) = 1;
        end
    end
end
legendre = [0 legendre(1:legendrelength-1)];   %生成legendre序列 

for k = 0 : WeilCodelength-1  
    WeilCode_p(k+1) = mod(sum([legendre(mod((k+intercept_pot_p-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff_p+intercept_pot_p-1), legendrelength)+1)]),2);
end
%将Weil码变换成极性码，0表示高电平‘+1’，1表示低电平‘-1’
for k = 1 : WeilCodelength
    if WeilCode_p(k) == 0
       WeilCode_p(k) = 1;
    else
       WeilCode_p(k) = -1;
    end    
end

end

