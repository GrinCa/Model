function FEmatrices = Modelv5_pattern(FEmatrices,listLHS,param,FILENAME)

ndof = size(FEmatrices.Nodes,1);
% regions
plate_region = 5;
ROOM_region = 6;
PML_region = 7;
%labels
in_plate_label = 2;
ext_plate_label = 3;
ROOMPML_label = 4;


Kr     = listLHS{1};
Ki     = listLHS{2}; % Stiffness matrix elastic domain
M      = listLHS{3}; % Mass matrix elastic domain
Hpmlr  = listLHS{4};
Hpmli  = listLHS{5};
Qpmlr  = listLHS{6};
Qpmli  = listLHS{7};
Ctmp1   = listLHS{8};
Ctmp2   = listLHS{9}; % coupling matrix

%connectivity
FEmatrices.connectivity = load(['Matrices/',FILENAME,'/connectivity_table.txt']);


% get the arrays of the nodes belonging to a given region/label
tab_region = get_regions([plate_region,ROOM_region,PML_region],ndof,FILENAME);
FEmatrices.plate_nodes  = find(tab_region(:,1));
FEmatrices.ROOM_nodes   = find(tab_region(:,2));
FEmatrices.PML_nodes    = find(tab_region(:,3));
FEmatrices.SIS          = find((tab_region(:,3)+tab_region(:,2)) >= 1);%semi-infinite-space

labels_cell = get_labels([in_plate_label,ext_plate_label,ROOMPML_label],...
                         FILENAME);
FEmatrices.PlateIn   = find(labels_cell{1});
FEmatrices.PlateExt  = find(labels_cell{2});
FEmatrices.ROOMPML   = find(labels_cell{3});


tab_plate = [];
for ii=1:length(FEmatrices.plate_nodes)
    for jj=1:3
        tab_plate = [tab_plate 3*(FEmatrices.plate_nodes(ii)-1)+jj];
    end
end




% plot3(FEmatrices.Nodes(FEmatrices.PlateExt,1),FEmatrices.Nodes(FEmatrices.PlateExt,2),FEmatrices.Nodes(FEmatrices.PlateExt,3),'+');

% plot3(FEmatrices.Nodes(FEmatrices.ROOM_nodes,1),FEmatrices.Nodes(FEmatrices.ROOM_nodes,2),FEmatrices.Nodes(FEmatrices.ROOM_nodes,3),'+');

% plot3(FEmatrices.Nodes(FEmatrices.PlateIn,1),FEmatrices.Nodes(FEmatrices.PlateIn,2),FEmatrices.Nodes(FEmatrices.PlateIn,3),'+');


for ii=1:length(FEmatrices.PlateExt)
    FEmatrices.indexPlateExt(ii) = length(tab_plate)+find(FEmatrices.SIS==FEmatrices.PlateExt(ii)); 
end

% indexing of the differents subspaces for partitionning
FEmatrices.indexu1        = 1:3:length(tab_plate);
FEmatrices.indexu2        = 2:3:length(tab_plate);
FEmatrices.indexu3        = 3:3:length(tab_plate);
FEmatrices.indexSIS       = (length(tab_plate)+1) : 1 : (length(tab_plate)+length(FEmatrices.SIS));

K = Kr(tab_plate,tab_plate) + 1i*Ki(tab_plate,tab_plate);
M = M(tab_plate,tab_plate);
Hpmlr = Hpmlr(FEmatrices.SIS,FEmatrices.SIS);
Hpmli = Hpmli(FEmatrices.SIS,FEmatrices.SIS);
Hpml = Hpmlr + 1i*Hpmli;
Qpmlr = Qpmlr(FEmatrices.SIS,FEmatrices.SIS);
Qpmli = Qpmli(FEmatrices.SIS,FEmatrices.SIS);
Qpml = Qpmlr + 1i*Qpmli;
C = sparse(size(M,1),size(Hpml,2));
C(FEmatrices.indexu3,:) = Ctmp2(FEmatrices.plate_nodes,FEmatrices.SIS);
FEmatrices.SurfIn_matrix = Ctmp1(FEmatrices.PlateIn, FEmatrices.PlateIn);
FEmatrices.SurfExt_matrix = Ctmp2(FEmatrices.PlateExt, FEmatrices.PlateExt);

Kglob = sparse([K -C;...
                sparse(size(Hpml,1),size(K,2)) Hpml]);
Mglob = sparse([M sparse(size(M,1),size(Qpml,2));...
                param.rho*C' Qpml/param.c0^2]);

FEmatrices.LHS = {Kglob,Mglob};

FEmatrices.size_system = size(Kglob,1);

end


function tab_region = get_regions(region_array,ndof,FILENAME)

connectivity_table = load(['Matrices/',FILENAME,'/connectivity_table.txt']);
region_element = load(['Matrices/',FILENAME,'/regions.txt']);

tab_region = zeros(ndof,length(region_array));

for ii=1:length(region_array)
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


