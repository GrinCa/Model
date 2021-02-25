%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Convergence of the solution                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this script is useful to test the convergence of the method. It uses the
%.geo file and change the size of the elements, in order to change the
%number of nodes.
%%%
% Note : this algorithm use version 3.0.6 of Gmsh and version 4.1 of
% FreeFem. The algorithm is likely to fail if the version of Gmsh is higher
% than 3.0.6

% Add folders for Mumps, WCAWE, and Mesh functions/files

Meshfold = 'Matrices';
WCAWEfold = 'WCAWE';
Mumps = 'Mumps';
DataMap = 'DataMap';
Derivatives = 'Derivatives';
STL = 'STL';
Config = 'Config';
Timing = 'Timing';
PrePost = 'PrePost';

addpath(genpath(strcat(pwd,'/',Meshfold)));
addpath(genpath(strcat(pwd,'/',WCAWEfold)));
addpath(genpath(strcat(pwd,'/',Mumps)));
addpath(genpath(strcat(pwd,'/',DataMap)));
addpath(genpath(strcat(pwd,'/',Derivatives)));
addpath(genpath(strcat(pwd,'/',STL)));
addpath(genpath(strcat(pwd,'/',Config)));
addpath(genpath(strcat(pwd,'/',Timing)));
addpath(genpath(strcat(pwd,'/',PrePost)));


%clear all;


mesh.file = 'Modelv5';

config = str2func([mesh.file '_config']);
[~, param] = config();


remesh = 1;

if remesh

% both arrays stand for size of elastic and acoustic nodes respectively.
% For each colum of both arrays, the main script will be run. Therefore it
% is possible to put as much as values that we desire, but the size of each
% array must be the same. 
% NOTE : If the size of elastic and acoustic is different, the mesh
% elements won't be regular, which implies that the mean quadratic pressure
% is no longer a good indicator.
sizemesh = [0.04];
file_sizemesh = fopen('sizemesh.txt','wt');

n_mesh = length(sizemesh(1,:));

% data to change values in the .geo file before the compilation into .msh
keywords = {'sizemesh = '};
n_key = length(keywords);

% Import existing .geo file
path = ['Geometry/',mesh.file,'/'];
fid = fopen([path,mesh.file,'.geo'],'rt');
dataMesh = fread(fid);
fclose(fid);
dataMesh = char(dataMesh.');

for kk=1:n_mesh
    update_files(param,mesh);
    %save sizemesh to acces in the main_multi.m file.
    fprintf(file_sizemesh,[num2str(sizemesh(1,kk)),'\n']);
    path_soft = getpath();
    % compile .geo file
    disp('Compile .geo file...');
    command = [path_soft.gmsh,' -3 ',path,mesh.file,'.geo ','-o ',path,mesh.file,'.msh',' -format msh2'];
    % version gmsh : 3.0.6, whiwh match with FreeFem++ 4.1. Antother
    % version could not work
    system(command);
    main_multi;
end
else 
main_multi;
end
fclose(file_sizemesh);





