function [FEmatrices,ndof,timing,flag] = get_matrices(timing,flag,mesh,matrix_names,param)

isupdated = 0; % isupdated=1 : recalculate freefrem++ script

for ii=1:length(matrix_names)
    matrix_names{ii} = strcat(mesh.file,'/',matrix_names{ii});
    if exist(matrix_names(ii),"file")
    else
        isupdated=1;
    end
end


%--------------------------------------------------------------------------
%Run FreeFem++ script, IF NEEDED
%--------------------------------------------------------------------------

if (isupdated||flag.rerun) % EDP updated and/or mesh updated
   t_0 = cputime;
   disp('************************');
   disp('*Rerun FreeFem++ script*');
   disp('************************');
   edpcommand = strcat('FreeFem++'," ",mesh.file,'.edp');
   system(edpcommand);
   timing.freefem = cputime-t_0;
   disp('*********************************************************');
   output = sprintf('[Get_matrices:infos] CPUtime for building of matrices %.4f s',timing.freefem);
   disp(output);
   disp('*********************************************************');
   flag.recalculated = 1;
end % if

%--------------------------------------------------------------------------
% Get matrices from the files
%--------------------------------------------------------------------------

listLHS = cell(1,length(matrix_names)); 

tic;
% Matrices of the FE problem
for ii=1:length(matrix_names)
    fid = fopen(matrix_names(ii),'rt');
    for jj=1:3
        line = fgets(fid);
        if jj==3
           values = str2num(strtrim(line));
           nrows = values(1);
           ncolums = values(2);
        end
    end
    fclose(fid);
    matrix_data = importdata(matrix_names(ii)," ",3);
    matrix_data = matrix_data.data;
    %FEmatrices.listLHS{ii} = sparse(matrix_data(:,1)+1,matrix_data(:,2)+1,matrix_data(:,3));
    listLHS{ii} = sparse([matrix_data(:,1);nrows-1]+1,[matrix_data(:,2);ncolums-1]+1,[matrix_data(:,3);0]);
end

% Nodes
Nodes = load(strcat(mesh.file,'/',"Nodes.txt"));
ndof = size(Nodes,1);
toc;

% RHS
RHSdata = importdata(strcat(mesh.file,'/',"RHS.txt")," ",3);
RHSdata = RHSdata.data;
FEmatrices.RHS = diag(sparse(RHSdata(:,1)+1,RHSdata(:,2)+1,RHSdata(:,3)));

%--------------------------------------------------------------------------
% return
%--------------------------------------------------------------------------

FEmatrices.Nodes = Nodes; 
FEmatrices = build_global(FEmatrices,listLHS,param,mesh.file);
FEmatrices.RHS(FEmatrices.size_system) = 0;
end



function FEmatrices = build_global(FEmatrices,listLHS,param,FILENAME)

ndof = size(FEmatrices.Nodes,1);

K = listLHS{1}; % Stiffness matrix elastic domain
M = listLHS{2}; % Mass matrix elastic domain
H = listLHS{3}; % Stiffness matrix acoustic domain
Q = listLHS{4}; % Mass matrix acoustic domain
Ctmp = listLHS{5}; % coupling matrix


% label of the different region of the mesh
region_labels = load(['Matrices/',FILENAME,'/labels.txt']); 
cavity_label = 1;
plate_label = 2;

% initialisation arrays of respective nodes
cavity_nodes = zeros(ndof,1);
plate_nodes = zeros(ndof,1);

for ii=1:length(region_labels)
    if region_labels(ii,2) == cavity_label
        cavity_nodes(region_labels(ii,1)+1) = 1;
    elseif region_labels(ii,2) == plate_label
        plate_nodes(region_labels(ii,1)+1) = 1;
    end
end

cavity_nodes = find(cavity_nodes);
plate_nodes = find(plate_nodes);



tab_plate = [];
tab_cavity = cavity_nodes;

for ii=1:length(plate_nodes)
    for jj=1:3
        tab_plate = [tab_plate 3*(plate_nodes(ii)-1)+jj];
    end
end

% save arrays of Nodes, needed for partionning
FEmatrices.cavity_nodes = cavity_nodes;
FEmatrices.plate_nodes = plate_nodes;
FEmatrices.field = find(FEmatrices.Nodes(:,1)<(1e-10));

% plot3(FEmatrices.Nodes(cavity_nodes,1),FEmatrices.Nodes(cavity_nodes,2),FEmatrices.Nodes(cavity_nodes,3),'+');
% plot3(FEmatrices.Nodes(plate_nodes,1),FEmatrices.Nodes(plate_nodes,2),FEmatrices.Nodes(plate_nodes,3),'+');
% plot3(FEmatrices.Nodes(FEmatrices.field,1),FEmatrices.Nodes(FEmatrices.field,2),FEmatrices.Nodes(FEmatrices.field,3),'+');

% indexing of the differents subspaces for partitionning
FEmatrices.indexu1 = 1:3:length(tab_plate);
FEmatrices.indexu2 = 2:3:length(tab_plate);
FEmatrices.indexu3 = 3:3:length(tab_plate);
FEmatrices.indexp = length(tab_plate)+1:1:(length(tab_plate)+length(tab_cavity));

FEmatrices.indexfield = 3*find(FEmatrices.Nodes(plate_nodes,1)<(1e-10));


K = K(tab_plate,tab_plate);
M = M(tab_plate,tab_plate);
H = H(tab_cavity,tab_cavity);
Q = Q(tab_cavity,tab_cavity);
C = sparse(size(M,1),size(H,2));
C(FEmatrices.indexu1,:) = Ctmp(plate_nodes,cavity_nodes);

Kglob = sparse([K -C;sparse(size(H,1),size(K,2)) H]);
Mglob = sparse([M sparse(size(C,1),size(C,2));param.rho*C' Q/param.c0^2]);

FEmatrices.LHS = {Kglob,Mglob};

% size of the "reduced" system < 4*ndof
FEmatrices.size_system = size(FEmatrices.LHS{1},1);

end

































