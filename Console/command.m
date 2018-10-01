path = '..\Datasets\Moghadam\';
files = dir([path '*.jpg']);
len = length(files);
for i= 1:len 
    old_name = files(i).name;           %���������ļ���  
    pos = strfind(old_name, 'myimage0');   %�ҵ�spec-��λ��  
    new_name = old_name(pos+8:end);       %��spec-��ʼֱ��������е��ַ�Ϊ�µ��ļ���  
    old_path = strcat(path,old_name);   %��ȡԭ���ļ�����·��  
    new_path = strcat(path,new_name);   %��ȡ���ļ�����·��  
    status = movefile(old_path,new_path); %ͨ���ƶ��ļ����������ļ�  
end