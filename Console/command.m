path = '..\Datasets\Moghadam\';
files = dir([path '*.jpg']);
len = length(files);
for i= 1:len 
    old_name = files(i).name;           %仅仅包括文件名  
    pos = strfind(old_name, 'myimage0');   %找到spec-的位置  
    new_name = old_name(pos+8:end);       %从spec-开始直到最后所有的字符为新的文件名  
    old_path = strcat(path,old_name);   %获取原来文件名及路径  
    new_path = strcat(path,new_name);   %获取新文件名及路径  
    status = movefile(old_path,new_path); %通过移动文件，重命名文件  
end