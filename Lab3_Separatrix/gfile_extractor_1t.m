function [flux,RBDRY,ZBDRY,NBDRY,R,Z,time,rdim,zdim]=gfile_extractor_1t(shot_number_test,start_efit_time,MagMesh)

%функци€ извлечени€ данных из Gfailов 

n=1;

% while(1)
    filename=strcat('gfiles\',num2str(shot_number_test),'\g0', num2str(shot_number_test),'.00',num2str(start_efit_time)); % им€ гфайла в формате g0shot_number_test.00time   
    fid = fopen(filename);
    if (fid==-1) return; end
    fprintf(1, strcat('reading file ',' g0', num2str(shot_number_test),'.00',num2str(start_efit_time),'\n'));

% EFITD    07/24/96      # 30095 , 121ms           3  65  65
%  .890000000E+00  .140000000E+01  .850000000E+00  .100000000E-01  .000000000E+00
%  .411241149E+00  .119873626E-01 -.121669609E-01 -.881654604E-02  .393396000E+00

    
    
    scan_cell_all=textscan(fid,'%s'); %извлечение данных из файла
    data=scan_cell_all{1,1}; % дата выстрела
    shot_number=str2num(data{4}); % номер выстрела
    time(n)=sscanf(data{6},'%fms'); %момент времени выстрела - дл€ каждого файла свой

    rdim=str2num(data{10}); % размер сетки по радиусу в метрах
    zdim=str2num(data{11}); % размер сетки по Z в метрах
    zmid=str2num(data{14}); % середина измерени€ по Z

   
    delay=15+55*5; %57;%57+65*65/5; %%   % ; %O_o в строке п€ть значений. сдвигаем на п€тьдес€т п€ть строк отнсительно п€тнадцатой
    %(зачем тут эти строки- одному дъ€волу известно)

    for i=1:MagMesh
        for j=1:MagMesh
            flux(i,j)=str2num(data{delay}); % после этого сдвига читаем поток это сетка 65х65 в каждой €чейке которой значение потока
            delay=delay+1;
        end
    end

    for i=1:length(data)
        if (strcmp(data{i},'NBDRY')) % количество точе сепаратриссы
            NBDRY(n)=str2num(data{i+2});
            break;
        end
    end

    for i=1:length(data)
        if (strcmp(data{i},'RBDRY')) %координаты сепаратриссы по радиусу
            i=i+1;
            for k=1:NBDRY(n) 
                RBDRY(n,k)=str2num(data{i+k});
            end
            break;
        end
    end

    for i=1:length(data)
        if (strcmp(data{i},'ZBDRY')) %координаты сепаратриссы по Z
            i=i+1;
            for k=1:NBDRY(n) 
                ZBDRY(n,k)=str2num(data{i+k});
            end
            break;
        end
    end

    fclose(fid);
  
%дл€ считывани€ по t
% start_efit_time=start_efit_time+1;
% n=n+1;    
% end

% неправильный расчет координатной сетки!
Z=0.5*zdim*(-MagMesh+1:2:MagMesh-1)/MagMesh; % расчитываем координатную сетку в метрах
R=0.5*rdim*(1:2:2*MagMesh-1)/MagMesh;

% rmin=str2num(data{13});
% zmin=str2num(data{11})/2;
% dz=
% for i=1:MagMesh
%     for j=1:MagMesh
%         Rr=rmin+i*dr;
%         Zz=zmin+i*dz;
%     end
% end

end

