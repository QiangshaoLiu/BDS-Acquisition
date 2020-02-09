function [CAcode,CAcode_delay] = CAgen_fs_delay(PRN,fs,delay,sim_t)
% generateCAcode.m generates one of the 32 GPS satellite C/A codes.
%
%版本：2
%修改内容：在第一版能产生一周期1023长度的CA码的基础上增加了在fs采样频率下sim_t ms的实际采样CA码样点值及延时后的采样值
%
% [CAcode_delay] = CAgen_new(PRN,fs,delay)
%
%   Inputs:
%       PRN             - PRN number of the sequence.
%       fs              -采样率
%       delay           -信号达到接收端时延时的码片数
%
%   Outputs:
%       CAcode          - a vector containing the desired C/A code sequence 
%                         (chips).
%       CAcode_delay    -fs采样频率下1ms对应的CA码延时delay码片后的采样值

%--- Make the code shift array. The shift depends on the PRN number -------
% The g2s vector holds the appropriate shift of the g2 code to generate
% the C/A code (ex. for SV#19 - use a G2 shift of g2s(19) = 471)
g2s = [  5,   6,   7,   8,  17,  18, 139, 140, 141, 251, ...
       252, 254, 255, 256, 257, 258, 469, 470, 471, 472, ...
       473, 474, 509, 512, 513, 514, 515, 516, 859, 860, ...
       861, 862];

%--- Pick right shift for the given PRN number ----------------------------
g2shift = g2s(PRN);

%--- Generate G1 code -----------------------------------------------------

%--- Initialize g1 output to speed up the function ---
g1 = zeros(1, 1023);
%--- Load shift register ---
reg = -1*ones(1, 10);

%--- Generate all G1 signal chips based on the G1 feedback polynomial -----
for i=1:1023
    g1(i)       = reg(10);
    saveBit     = reg(3)*reg(10);
    reg(2:10)   = reg(1:9);
    reg(1)      = saveBit;
end

%--- Generate G2 code -----------------------------------------------------

%--- Initialize g2 output to speed up the function ---
g2 = zeros(1, 1023);
%--- Load shift register ---
reg = -1*ones(1, 10);

%--- Generate all G2 signal chips based on the G2 feedback polynomial -----
for i=1:1023
    g2(i)       = reg(10);
    saveBit     = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
    reg(2:10)   = reg(1:9);
    reg(1)      = saveBit;
end

%--- Shift G2 code --------------------------------------------------------
%The idea: g2 = concatenate[ g2_right_part, g2_left_part ];
g2 = [g2(1023-g2shift+1 : 1023), g2(1 : 1023-g2shift)];

%--- Form single sample C/A code by multiplying G1 and G2 -----------------
CAcode = -(g1 .* g2);


t_delay=delay/1023:1e3/fs:sim_t-1e3/fs+delay/1023;  %延迟delay的码片的一个码周期的对应时间，单位ms;零延迟的话对应delay=0
l=length(t_delay);
n_delay=fix(mod(t_delay,1)*1023)+1;             %延迟delay码片每次采样点落在第n_delay码片上
CAcode_delay=zeros(1,l);
for i=1:l
    CAcode_delay(i)=CAcode(n_delay(i));
end


