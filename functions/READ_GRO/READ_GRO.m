% Program by Sumith YD to READ gro file of water

function [gro_tensor, box_dimensions] = READ_GRO(filename)
%[filename, pathname] = uigetfile({'*.gro'},'File Selector');
% gro_tensor has five colums: molecule nr, atom nr, x, y, z
% box_dimensions corresponds to the last line of the .gro file

fid = fopen(filename);
tline = fgetl(fid);
line_NUM=1;
while ischar(tline)             % find total number of lines of the .gro file for initialization
    tline = fgetl(fid);
    line_NUM = line_NUM +1;
end
TOTAL_LINES = line_NUM-1;
fclose(fid);
napht = zeros(line_NUM-4,5);
box_dim = zeros(9,1);
fid = fopen(filename);
tline = fgetl(fid);             % reads line number fid
line_NUM=1;
naphtCount = 1;
while ischar(tline)
    t = strsplit(tline, ' ');
    if(line_NUM > 2 && line_NUM < TOTAL_LINES)
        m = char(t(2));
        %moleculeType=cellstr(m(length(m)-2:length(m)));
        %if(strcmp(moleculeType(1),'NAP'))
            napht(naphtCount,1) = str2double(m(1:length(m)-3));
            napht(naphtCount,2) = str2double(t(4));
            napht(naphtCount,3) = str2double(t(5));
            napht(naphtCount,4) = str2double(t(6));
            napht(naphtCount,5) = str2double(t(7));
            naphtCount = naphtCount + 1;
        %end
    end
    if (line_NUM == TOTAL_LINES)
        if length(t)==3 || length(t)==9
            for i=1:9
                if i<=length(t)
                    box_dim(i) = str2double(t(i));
                else
                    box_dim(i) = 0;
                end
            end
        end
        if length(t)==4
            for i=1:9
                if i<=3
                    box_dim(i) = str2double(t(i+1));
                else
                    box_dim(i) = 0;
                end
            end
        end
        if length(t)==10
            for i=1:9
                box_dim(i) = str2double(t(i+1));
            end
        end
    end
    tline = fgetl(fid);
    line_NUM = line_NUM +1;
end
fclose(fid);
gro_tensor = napht;
box_dimensions = box_dim;
%msgbox('File read into COORDINATES');
end