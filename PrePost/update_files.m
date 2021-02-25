function update_files(param,mesh)

% this function modifies the content of the .edp and .geo files. It takes
% the struct param with fields EDP and GEO that contain the keywords 
% and values fields. These fields have to be defined in the Config file.


% Change EDP file
EDP_file = ['EDP/',mesh.file,'/',mesh.file,'.edp'];
EDP_fid = fopen(EDP_file,'r+');
updated_text = [];
while ~feof(EDP_fid)
    line = fgets(EDP_fid);
    for ii=1:size(param.EDP.keywords,1)
        idx = strfind(param.EDP.keywords{ii,1},'=');
        if strncmpi(line,param.EDP.keywords{ii,1},idx(1))
            line = sprintf(param.EDP.keywords{ii,1},param.EDP.keywords{ii,2});
            break;
        end
    end
    updated_text = [updated_text line];
end
fclose(EDP_fid);
%save modified content
fid = fopen(EDP_file,'w+');
fprintf(fid,updated_text);
fclose(fid);


% Change GEO file
GEO_file = ['Geometry/',mesh.file,'/',mesh.file,'.geo'];
GEO_fid = fopen(GEO_file,'r+');
updated_text = [];
while ~feof(GEO_fid)
    line = fgets(GEO_fid);
    for ii=1:size(param.GEO.keywords,1)
        idx = strfind(param.GEO.keywords{ii,1},'=');
        if strncmpi(line,param.GEO.keywords{ii,1},idx(1))
            line = sprintf(param.GEO.keywords{ii,1},param.GEO.keywords{ii,2});
        end
    end
    updated_text = [updated_text line];
end
fclose(GEO_fid);
%save modified content
fid = fopen(GEO_file,'w+');
fprintf(fid,updated_text);
fclose(fid);

end
