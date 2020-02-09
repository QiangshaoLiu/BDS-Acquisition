function [ WeilCode_d ] = generatedatacode( PRN )
%GENERATECODE Summary of this function goes here
%   Detailed explanation goes here
%产生PRN号为1至10的主码序列

switch PRN
    case (1)
        phase_diff_d = 2678;
        intercept_pot_d = 699;
    case (2)
        phase_diff_d = 4802;
        intercept_pot_d = 694;
    case (3)
        phase_diff_d = 958;
        intercept_pot_d = 7318;
    case (4)
        phase_diff_d = 859;
        intercept_pot_d = 2127;
    case (5)
        phase_diff_d = 3843;
        intercept_pot_d = 715;
    case (6)
        phase_diff_d = 2232;
        intercept_pot_d = 6682;
    case (7)
        phase_diff_d = 124;
        intercept_pot_d = 7850;
    case (8)
        phase_diff_d = 4352;
        intercept_pot_d = 5495;
    case(9)
        phase_diff_d = 1816;
        intercept_pot_d = 1162;
    case (10)
        phase_diff_d = 1126;
        intercept_pot_d = 7682;
    case (11)
        phase_diff_d = 1860;
        intercept_pot_d = 6792;
    case (12)
        phase_diff_d = 4800;
        intercept_pot_d = 9973;    
    case (13)
        phase_diff_d = 2267;
        intercept_pot_d = 6596;     
    case (14)
        phase_diff_d = 424;
        intercept_pot_d = 2092;    
    case (15)
        phase_diff_d = 4192;
        intercept_pot_d = 19;    
    case (16)
        phase_diff_d = 4333;
        intercept_pot_d = 10151;     
    case (17)
        phase_diff_d = 2656;
        intercept_pot_d = 6297;
    case (18)
        phase_diff_d = 4148;
        intercept_pot_d = 5766;
    case (19)
        phase_diff_d = 243;
        intercept_pot_d = 2359;
    case (20)
        phase_diff_d = 1330;
        intercept_pot_d = 7136;
    case (21)
        phase_diff_d = 1593;
        intercept_pot_d = 1706;    
    case (22)
        phase_diff_d = 1470;
        intercept_pot_d = 2128;    
    case (23)
        phase_diff_d = 882;
        intercept_pot_d = 6827;    
    case (24)
        phase_diff_d = 3202;
        intercept_pot_d = 693;    
    case (25)
        phase_diff_d = 5095;
        intercept_pot_d = 9729;    
    case (26)
        phase_diff_d = 2546;
        intercept_pot_d = 1620;    
    case (27)
        phase_diff_d = 1733;
        intercept_pot_d = 6805;    
    case (28)
        phase_diff_d = 4795;
        intercept_pot_d = 534;    
    case (29)
        phase_diff_d = 4577;
        intercept_pot_d = 712; 
    case (30)
        phase_diff_d = 1627;
        intercept_pot_d = 1929;
    case (31)
        phase_diff_d = 3638;
        intercept_pot_d = 5355;    
    case (32)
        phase_diff_d = 2553;
        intercept_pot_d = 6139;    
    case (33)
        phase_diff_d = 3646;
        intercept_pot_d = 6339;    
    case (34)
        phase_diff_d = 1087;
        intercept_pot_d = 1470;     
    case (35)
        phase_diff_d = 1843;
        intercept_pot_d = 6867; 
    case (36)
        phase_diff_d = 216;
        intercept_pot_d = 7851;     
    case (37)
        phase_diff_d = 2245;
        intercept_pot_d = 1162;    
    case (38)
        phase_diff_d = 726;
        intercept_pot_d = 7659; 
    case (39)
        phase_diff_d = 1966;
        intercept_pot_d = 1156;    
    case (40)
        phase_diff_d = 670;
        intercept_pot_d = 2672;    
    case (41)
        phase_diff_d = 4130;
        intercept_pot_d = 6043;    
    case (42)
        phase_diff_d = 53;
        intercept_pot_d = 2862;    
    case (43)
        phase_diff_d = 4830;
        intercept_pot_d = 180;    
    case (44)
        phase_diff_d = 182;
        intercept_pot_d = 2663;    
    case (45)
        phase_diff_d = 2181;
        intercept_pot_d = 6940;    
    case (46)
        phase_diff_d = 2006;
        intercept_pot_d = 1645;
    case (47)
        phase_diff_d = 1080;
        intercept_pot_d = 1582;    
    case (48)
        phase_diff_d = 2288;
        intercept_pot_d = 951;    
    case (49)
        phase_diff_d = 2027;
        intercept_pot_d = 6878;     
    case (50)
        phase_diff_d = 271;
        intercept_pot_d = 7701;    
    case (51)
        phase_diff_d = 915;
        intercept_pot_d = 1823; 
    case (52)
        phase_diff_d = 497;
        intercept_pot_d = 2391;    
    case (53)
        phase_diff_d = 139;
        intercept_pot_d = 2606;    
    case (54)
        phase_diff_d = 3693;
        intercept_pot_d = 822;    
    case (55)
        phase_diff_d = 2054;
        intercept_pot_d = 6403;    
    case (56)
        phase_diff_d = 4342;
        intercept_pot_d = 239;    
    case (57)
        phase_diff_d = 3342;
        intercept_pot_d = 442;   
    case (58)
        phase_diff_d = 2592;
        intercept_pot_d = 6769;  
    case (59)
        phase_diff_d = 1007;
        intercept_pot_d = 2560;
    case (60)
        phase_diff_d = 310;
        intercept_pot_d = 2502;    
    case (61)
        phase_diff_d = 4203;
        intercept_pot_d = 5072;    
    case (62)
        phase_diff_d = 455;
        intercept_pot_d = 7268;    
    case (63)
        phase_diff_d = 4318;
        intercept_pot_d = 341;             
    otherwise
        phase_diff_d = 2678;
        intercept_pot_d = 699;    %超出范围，默认PRN=1
end 

legendrelength = 10243;
WeilCodelength = 10230;
legendre = zeros(1, legendrelength);
WeilCode_d = zeros(1, WeilCodelength);

for k = 1 : legendrelength-1 
    for x = 1 : (legendrelength-1)/2
        if mod(k,legendrelength) == mod(x^2, legendrelength)
            legendre(k) = 1;
        end
    end
end
legendre = [0 legendre(1:legendrelength-1)];   %生成legendre序列 

for k = 0 : WeilCodelength-1  
    WeilCode_d(k+1) = mod(sum([legendre(mod((k+intercept_pot_d-1), legendrelength)+1),... 
        legendre(mod((k+phase_diff_d+intercept_pot_d-1), legendrelength)+1)]),2);
end
%将Weil码变换成极性码，0表示高电平‘+1’，1表示低电平‘-1’
for k = 1 : WeilCodelength
    if WeilCode_d(k) == 0
       WeilCode_d(k) = 1;
    else
       WeilCode_d(k) = -1;
    end    
end

end