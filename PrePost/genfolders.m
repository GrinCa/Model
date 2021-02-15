function param = genfolders(mesh,param)
Filename = mesh.file;

%--------------------------------------------------------------------------
% Matrices
%--------------------------------------------------------------------------

path1 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']',...
         '[' num2str(int16(180/pi*param.theta_range(1))) '_' num2str(int16(180/pi*param.theta_range(2))) ']'];

if exist(['Matrices/',Filename],'dir') == 0
    path = ['Matrices/',Filename];
    mkdir(path);
end

if exist(['Matrices/',Filename,'/',param.path1],'dir') == 0
    path = ['Matrices/',Filename,'/',param.path1];
    mkdir(path);
end

if exist(['Mesh/',Filename],'dir') == 0
    path = ['Mesh/',Filename];
    mkdir(path);
end

if exist(['EDP/',Filename],'dir') == 0
    path = ['EDP/',Filename];
    mkdir(path);
end


for ii=1:length(param.vecfreqrange)
    for jj=1:length(param.vecthetarange)
        if exist(['Matrices/',Filename,'/',param.path2],'dir') == 0
            path = ['Matrices/',Filename,'/',param.path2];
            mkdir(path);
        end
    end
end

%--------------------------------------------------------------------------
% Geometry : copy .msh from Geometry/filename into the main folder
%--------------------------------------------------------------------------

try
    command = ['cp ','Geometry/',Filename,'/',Filename,'.msh',' ','Mesh/',Filename,'/',Filename,'.msh'];
    system(command);
catch
    disp(['[genfolders] No .msh file in the folder : /Geometry/',Filename]);
    return;
end
%--------------------------------------------------------------------------
% DataMap
%--------------------------------------------------------------------------

if exist(['DataMap/',Filename]) == 0
    command = ['DataMap/',Filename];
    system(['mkdir ' command]);
end

param.interval_detail_str = '';

for n=1:2 % frequency + angle
    for ii=1:length(param.interval_construct{n})
        if length(param.interval_construct{n}{ii}) == 1 % loop on frequency/angle sub interval
            param.interval_detail_str = strcat(param.interval_detail_str,['[' num2str(param.interval_construct{n}{ii}) ']']);
        else
            for kk=1:length(param.interval_construct{n}{ii})
                if kk==1
                    param.interval_detail_str = strcat(param.interval_detail_str,['[' num2str(param.interval_construct{n}{ii}(kk))]);
                elseif kk==length(param.interval_construct{1}{ii})
                    param.interval_detail_str = strcat(param.interval_detail_str,[num2str(param.interval_construct{n}{ii}(kk)) ']']);
                    break;
                else
                    param.interval_detail_str = strcat(param.interval_detail_str,num2str(param.interval_construct{n}{ii}(kk)));
                end
                param.interval_detail_str = strcat(param.interval_detail_str,'_');
            end
        end
    end
    if n==1
        param.interval_detail_str = strcat(param.interval_detail_str,'#');
    end
end


    

end