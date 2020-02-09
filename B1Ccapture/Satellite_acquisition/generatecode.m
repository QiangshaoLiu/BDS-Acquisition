function [ WeilCode_p ] = generatecode( PRN )
%GENERATECODE Summary of this function goes here
%   Detailed explanation goes here
%产生PRN号为1至63的导频分量伪码主码序列

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
    case (11)
        phase_diff_p = 2715;
        intercept_pot_p = 1402;
    case (12)
        phase_diff_p = 4408;
        intercept_pot_p = 5557;
    case (13)
        phase_diff_p = 3160;
        intercept_pot_p = 5764;
    case (14)
        phase_diff_p = 2796;
        intercept_pot_p = 1073;
    case (15)
        phase_diff_p = 459;
        intercept_pot_p = 7001;
    case (16)
        phase_diff_p = 3594;
        intercept_pot_p = 5910;
    case (17)
        phase_diff_p = 4813;
        intercept_pot_p = 10060;
    case (18)
        phase_diff_p = 586;
        intercept_pot_p = 2710;
    case (19)
        phase_diff_p = 1428;
        intercept_pot_p = 1546;
    case (20)
        phase_diff_p = 2371;
        intercept_pot_p = 6887;
    case (21)
        phase_diff_p = 2285;
        intercept_pot_p = 1883;
    case (22)
        phase_diff_p = 3377;
        intercept_pot_p = 5613;
    case (23)
        phase_diff_p = 4965;
        intercept_pot_p = 5062;
    case (24)
        phase_diff_p = 3779;
        intercept_pot_p = 1038;
    case (25)
        phase_diff_p = 4547;
        intercept_pot_p = 10170;
    case (26)
        phase_diff_p = 1646;
        intercept_pot_p = 6484;
    case (27)
        phase_diff_p = 1430;
        intercept_pot_p = 1718;
    case (28)
        phase_diff_p = 607;
        intercept_pot_p = 2535;
    case (29)
        phase_diff_p = 2118;
        intercept_pot_p = 1158;
    case (30)
        phase_diff_p = 4709;
        intercept_pot_p = 526;
    case (31)
        phase_diff_p = 1149;
        intercept_pot_p = 7331;
    case (32)
        phase_diff_p = 3283;
        intercept_pot_p = 5844;
    case (33)
        phase_diff_p = 2473;
        intercept_pot_p = 6423;
    case (34)
        phase_diff_p = 1006;
        intercept_pot_p = 6968;
    case (35)
        phase_diff_p = 3670;
        intercept_pot_p = 1280;
    case (36)
        phase_diff_p = 1817;
        intercept_pot_p = 1838;
    case (37)
        phase_diff_p = 771;
        intercept_pot_p = 1989;
    case (38)
        phase_diff_p = 2173;
        intercept_pot_p = 6468;
    case (39)
        phase_diff_p = 740;
        intercept_pot_p = 2091;
    case (40)
        phase_diff_p = 1433;
        intercept_pot_p = 1581;
    case (41)
        phase_diff_p = 2458;
        intercept_pot_p = 1453;
    case (42)
        phase_diff_p = 3459;
        intercept_pot_p = 6252;
    case (43)
        phase_diff_p = 2155;
        intercept_pot_p = 7122;
    case (44)
        phase_diff_p = 1205;
        intercept_pot_p = 7711;
    case (45)
        phase_diff_p = 413;
        intercept_pot_p = 7216;
    case (46)
        phase_diff_p = 874;
        intercept_pot_p = 2113;
    case (47)
        phase_diff_p = 2463;
        intercept_pot_p = 1095;
    case (48)
        phase_diff_p = 1106;
        intercept_pot_p = 1628;
    case (49)
        phase_diff_p = 1590;
        intercept_pot_p = 1713;
    case (50)
        phase_diff_p = 3873;
        intercept_pot_p = 6102;
    case (51)
        phase_diff_p = 4026;
        intercept_pot_p = 6123;
    case (52)
        phase_diff_p = 4272;
        intercept_pot_p = 6070;
    case (53)
        phase_diff_p = 3556;
        intercept_pot_p = 1115;
    case (54)
        phase_diff_p = 128;
        intercept_pot_p = 8047;
    case (55)
        phase_diff_p = 1200;
        intercept_pot_p = 6795;
    case (56)
        phase_diff_p = 130;
        intercept_pot_p = 2575;
    case (57)
        phase_diff_p = 4494;
        intercept_pot_p = 53;
    case (58)
        phase_diff_p = 1871;
        intercept_pot_p = 1729;
    case (59)
        phase_diff_p = 3073;
        intercept_pot_p = 6388;
    case (60)
        phase_diff_p = 4386;
        intercept_pot_p = 682;
    case (61)
        phase_diff_p = 4098;
        intercept_pot_p = 5565;
    case (62)
        phase_diff_p = 1923;
        intercept_pot_p = 7160;
    case (63)
        phase_diff_p = 1176;
        intercept_pot_p = 2277;
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