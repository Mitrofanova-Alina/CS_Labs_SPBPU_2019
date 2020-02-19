function [flux,RBDRY,ZBDRY,NBDRY,R,Z,time,rdim,zdim]=gfile_extractor_1t(shot_number_test,start_efit_time,MagMesh)

%������� ���������� ������ �� Gfail�� 

n=1;

% while(1)
    filename=strcat('gfiles/',num2str(shot_number_test),'/g0', num2str(shot_number_test),'.00',num2str(start_efit_time)); % ��� ������ � ������� g0shot_number_test.00time   
    fid = fopen(filename);
    if (fid==-1) return; end
    fprintf(1, strcat('reading file ',' g0', num2str(shot_number_test),'.00',num2str(start_efit_time),'\n'));

% EFITD    07/24/96      # 30095 , 121ms           3  65  65
%  .890000000E+00  .140000000E+01  .850000000E+00  .100000000E-01  .000000000E+00
%  .411241149E+00  .119873626E-01 -.121669609E-01 -.881654604E-02  .393396000E+00

    
    
    scan_cell_all=textscan(fid,'%s'); %���������� ������ �� �����
    data=scan_cell_all{1,1}; % ���� ��������
    shot_number=str2num(data{4}); % ����� ��������
    time(n)=sscanf(data{6},'%fms'); %������ ������� �������� - ��� ������� ����� ����

    rdim=str2num(data{10}); % ������ ����� �� ������� � ������
    zdim=str2num(data{11}); % ������ ����� �� Z � ������
    zmid=str2num(data{14}); % �������� ��������� �� Z

   
    delay=15+55*5; %57;%57+65*65/5; %%   % ; %O_o � ������ ���� ��������. �������� �� ��������� ���� ����� ����������� �����������
    %(����� ��� ��� ������- ������ ������� ��������)

    for i=1:MagMesh
        for j=1:MagMesh
            flux(i,j)=str2num(data{delay}); % ����� ����� ������ ������ ����� ��� ����� 65�65 � ������ ������ ������� �������� ������
            delay=delay+1;
        end
    end

    for i=1:length(data)
        if (strcmp(data{i},'NBDRY')) % ���������� ���� ������������
            NBDRY(n)=str2num(data{i+2});
            break;
        end
    end

    for i=1:length(data)
        if (strcmp(data{i},'RBDRY')) %���������� ������������ �� �������
            i=i+1;
            for k=1:NBDRY(n) 
                RBDRY(n,k)=str2num(data{i+k});
            end
            break;
        end
    end

    for i=1:length(data)
        if (strcmp(data{i},'ZBDRY')) %���������� ������������ �� Z
            i=i+1;
            for k=1:NBDRY(n) 
                ZBDRY(n,k)=str2num(data{i+k});
            end
            break;
        end
    end

    fclose(fid);
  
%��� ���������� �� t
% start_efit_time=start_efit_time+1;
% n=n+1;    
% end

% ������������ ������ ������������ �����!
Z=0.5*zdim*(-MagMesh+1:2:MagMesh-1)/MagMesh; % ����������� ������������ ����� � ������
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
