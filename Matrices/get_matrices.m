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


% Matrices of the FE problem
tic;

for ii=1:length(matrix_names)
    nrows = 0;
    ncolums = 0;
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
    listLHS{ii} = sparse([matrix_data(:,1);nrows-1]+1,[matrix_data(:,2);ncolums-1]+1,[matrix_data(:,3);0]);
end
toc;
% Nodes
Nodes = load(strcat(mesh.file,'/',"Nodes.txt"));
ndof = size(Nodes,1);



% RHS
% RHSdata = importdata(strcat(mesh.file,'/',"RHS.txt")," ",3);
% RHSdata = RHSdata.data;
% FEmatrices.RHS = diag(sparse(RHSdata(:,1)+1,RHSdata(:,2)+1,RHSdata(:,3)));

%--------------------------------------------------------------------------
% return
%--------------------------------------------------------------------------

FEmatrices.Nodes = Nodes; 
FEmatrices = build_global(FEmatrices,listLHS,param,mesh.file);
FEmatrices = build_BGField(FEmatrices, param);

end


function FEmatrices = build_global(FEmatrices,listLHS,param,FILENAME)

FEmatrices = Model_pattern(FEmatrices,listLHS,param,FILENAME);

end


% this following function build the matrices of the system for the plate
% case only. For each new case, we have to build a new function, a pattern
% for partitionning the system.
function FEmatrices = Model_pattern(FEmatrices,listLHS,param,FILENAME)

ndof = size(FEmatrices.Nodes,1);
% regions
cavity_region = 1;
plate_region = 2;
background_region = 3;
PML_region = 4;
%labels
PlateCavity_label = 5;
PlateBG_label = 6;
BGPML_label = 8;


K     = listLHS{1}; % Stiffness matrix elastic domain
M     = listLHS{2}; % Mass matrix elastic domain
Hbg   = listLHS{3}; % Stiffness matrix acoustic domain
Qbg   = listLHS{4}; % Mass matrix acoustic domain
Hcav  = listLHS{5};
Qcav  = listLHS{6};
Hpmlr = listLHS{7};
Hpmli = listLHS{8};
Qpmlr = listLHS{9};
Qpmli = listLHS{10};
C1tmp = listLHS{11}; % coupling matrix
C2tmp = listLHS{12};


% get the arrays of the nodes belonging to a given region/label
tab_region = get_regions([cavity_region,plate_region,background_region,PML_region],ndof,FILENAME);
FEmatrices.cavity_nodes     = find(tab_region(:,1));
FEmatrices.plate_nodes      = find(tab_region(:,2));
FEmatrices.BG_nodes         = find(tab_region(:,3));
FEmatrices.PML_nodes        = find(tab_region(:,4));
FEmatrices.BG_PML_nodes     = find((tab_region(:,3) + tab_region(:,4)) >= 1);

labels_cell = get_labels([PlateCavity_label,PlateBG_label,BGPML_label],...
                         FILENAME);
FEmatrices.PlateCavity_nodes    = find(labels_cell{1});
FEmatrices.PlateBG_nodes        = find(labels_cell{2});
FEmatrices.BGPML_nodes          = find(labels_cell{3});

tab_plate = [];
for ii=1:length(FEmatrices.plate_nodes)
    for jj=1:3
        tab_plate = [tab_plate 3*(FEmatrices.plate_nodes(ii)-1)+jj];
    end
end


%plot3(FEmatrices.Nodes(cavity_nodes,1),FEmatrices.Nodes(cavity_nodes,2),FEmatrices.Nodes(cavity_nodes,3),'+');
% plot3(FEmatrices.Nodes(plate_nodes,1),FEmatrices.Nodes(plate_nodes,2),FEmatrices.Nodes(plate_nodes,3),'+');
% plot3(FEmatrices.Nodes(FEmatrices.field,1),FEmatrices.Nodes(FEmatrices.field,2),FEmatrices.Nodes(FEmatrices.field,3),'+');

% indexing of the differents subspaces for partitionning
FEmatrices.indexu1        = 1:3:length(tab_plate);
FEmatrices.indexu2        = 2:3:length(tab_plate);
FEmatrices.indexu3        = 3:3:length(tab_plate);
FEmatrices.indexP_BG_PML  = (length(tab_plate)+1) : 1 : (length(tab_plate)+length(FEmatrices.BG_PML_nodes));
FEmatrices.indexP_CAVITY  = (length(tab_plate)+length(FEmatrices.BG_PML_nodes)+1) : 1 : ...
                            (length(tab_plate)+length(FEmatrices.BG_PML_nodes)+length(FEmatrices.cavity_nodes));
FEmatrices.indexP_BG = zeros(length(FEmatrices.BG_nodes),1);
for ii=1:length(FEmatrices.BG_nodes)
    FEmatrices.indexP_BG(ii) = length(tab_plate)+find(FEmatrices.BG_PML_nodes==FEmatrices.BG_nodes(ii)); 
end

for ii=1:length(FEmatrices.BGPML_nodes)
    FEmatrices.indexP_BGPML(ii) = length(tab_plate)+find(FEmatrices.BG_PML_nodes==FEmatrices.BGPML_nodes(ii)); 
end


K = K(tab_plate,tab_plate);
M = M(tab_plate,tab_plate);
FEmatrices.Hbg = Hbg(FEmatrices.BG_nodes,FEmatrices.BG_nodes);%needed for the calculation of the RHS
FEmatrices.Qbg = Qbg(FEmatrices.BG_nodes,FEmatrices.BG_nodes);%needed for the calculation of the RHS
Hcav = Hcav(FEmatrices.cavity_nodes,FEmatrices.cavity_nodes);
Qcav = Qcav(FEmatrices.cavity_nodes,FEmatrices.cavity_nodes);
Hpmlr = Hpmlr(FEmatrices.BG_PML_nodes,FEmatrices.BG_PML_nodes);
Hpmli = Hpmli(FEmatrices.BG_PML_nodes,FEmatrices.BG_PML_nodes);
Hpml = Hpmlr + 1i*Hpmli;
Qpmlr = Qpmlr(FEmatrices.BG_PML_nodes,FEmatrices.BG_PML_nodes);
Qpmli = Qpmli(FEmatrices.BG_PML_nodes,FEmatrices.BG_PML_nodes);
Qpml = Qpmlr + 1i*Qpmli;
C1 = sparse(size(M,1),size(Hpml,2));
C2 = sparse(size(M,1),size(Hcav,2));
C1(FEmatrices.indexu1,:) = C1tmp(FEmatrices.plate_nodes,FEmatrices.BG_PML_nodes);
C2(FEmatrices.indexu1,:) = C2tmp(FEmatrices.plate_nodes,FEmatrices.cavity_nodes);
FEmatrices.C1 = C1tmp(FEmatrices.PlateBG_nodes, FEmatrices.PlateBG_nodes);
FEmatrices.C2 = C2tmp(FEmatrices.PlateCavity_nodes, FEmatrices.PlateCavity_nodes);

Kglob = sparse([K -C1 -C2;...
                sparse(size(Hpml,1),size(K,2)) Hpml sparse(size(Hpml,1),size(Hcav,2));...
                sparse(size(Hcav,1),size(K,2)+size(Hpml,2)) Hcav]);
Mglob = sparse([M sparse(size(M,1),size(Qpml,2)+size(Qcav,2));...
                param.rho*C1' Qpml/param.c0^2 sparse(size(Qpml,1),size(Qcav,2));...
                param.rho*C2' sparse(size(Qcav,1),size(Qpml,2)) Qcav/param.c0^2]);

FEmatrices.LHS = {Kglob,Mglob};

FEmatrices.size_system = size(Kglob,1);

% FEmatrices.indexfield is the index of the nodes (of the RHS vector)  where is applied the BG
% pressure field
FEmatrices.field = FEmatrices.BG_nodes;
FEmatrices.indexfield = FEmatrices.indexP_BG;
end

function FEmatrices = build_BGField(FEmatrices,param)
ndof = size(FEmatrices.Nodes,1);

% Cell of RHS
FEmatrices.RHS_BG = zeros(FEmatrices.size_system,param.nfreq,param.ntheta);

%BG_nodes = 1:1:FEmatrices.size_system;
Field_nodes = FEmatrices.field;

xbg = FEmatrices.Nodes(Field_nodes,1);
ybg = FEmatrices.Nodes(Field_nodes,2);

FEmatrices.BG_pressure = zeros(ndof,param.nfreq,param.ntheta);

P0 = 1;

for ii=1:param.nfreq
    for jj=1:param.ntheta
        k = 2*pi*param.freq(ii)/param.c0;
        BG_Pressure_tmp = P0*exp(1i*k*(xbg*cos(param.theta(jj))+ybg*sin(param.theta(jj))));%propagation (+x,+y), convention exp(-1i*k*x)
        FEmatrices.BG_pressure(Field_nodes,ii,jj) = BG_Pressure_tmp;
        Z = FEmatrices.Hbg - (2*pi*param.freq(ii)/param.c0)^2*FEmatrices.Qbg;
        U_inc = zeros(FEmatrices.size_system,1);
        U_inc(FEmatrices.indexfield) = -Z*BG_Pressure_tmp;
        U_inc(FEmatrices.indexP_BGPML) = 0;
        FEmatrices.RHS_BG(:,ii,jj) = U_inc;
    end
end
end


function tab_region = get_regions(region_array,ndof,FILENAME)

connectivity_table = load(['Matrices/',FILENAME,'/connectivity_table.txt']);
region_element = load(['Matrices/',FILENAME,'/regions.txt']);

tab_region = zeros(ndof,length(region_array));

for ii=region_array
    id_elements = find(region_element == region_array(ii));
    for jj=1:length(id_elements)
        tab_region(connectivity_table(id_elements(jj),:)+1,ii) = 1;
    end
end

end

function labels_cell = get_labels(label_number,FILENAME)

tab_labels = load(['Matrices/',FILENAME,'/labels.txt']);
labels_cell = cell(1,length(label_number));

for ii=1:length(label_number)
    jj = find(tab_labels(1,:) == label_number(ii));
    labels_cell{ii} = tab_labels(2:end,jj);
end

end





























