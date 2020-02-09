% 本程序用于生成正弦波查表所需的数据
% 生成 ROM 的 .coe文件
fid = fopen('cos_rom.txt','w');

fprintf(fid,'MEMORY_INITIALIZATION_RADIX = 10;\n');

fprintf(fid,'MEMORY_INITIALIZATION_VECTOR =\n');

for i = 0:1:pi/2*100

    y = cos(i/100);

    rom =floor( y * 2^12);

    if i == 157

        fprintf(fid,'%d;',rom);

    else

        fprintf(fid,'%d,',rom);

    end

   

    if mod(i,10)==0 && i ~= 0

        fprintf(fid,'\n');

    end

end

fclose(fid);