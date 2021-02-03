function FEmatrices = Modelv3_pattern(FEmatrices,listLHS,param,FILENAME)

ndof = size(FEmatrices.Nodes,1);
% regions
plate_region = 5;
BG_region = 6;
PML_region = 7;
%labels
ext_plate_label = 2;
in_plate_label = 3;
BGPML_label = 4;


Kr     = listLHS{1};
Ki     = listLHS{2}; % Stiffness matrix elastic domain
M      = listLHS{3}; % Mass matrix elastic domain
Hbg    = listLHS{4}; % Stiffness matrix acoustic domain
Qbg    = listLHS{5}; % Mass matrix acoustic domain
Hpmlr  = listLHS{6};
Hpmli  = listLHS{7};
Qpmlr  = listLHS{8};
Qpmli  = listLHS{9};
Ctmp   = listLHS{10}; % coupling matrix


%connectivity
FEmatrices.connectivity = load(['Matrices/',FILENAME,'/connectivity_table.txt']);


% get the arrays of the nodes belonging to a given region/label
tab_region = get_regions([plate_region,BG_region,PML_region],ndof,FILENAME);
FEmatrices.plate_nodes  = find(tab_region(:,1));
FEmatrices.BG_nodes     = find(tab_region(:,2));
FEmatrices.PML_nodes    = find(tab_region(:,3));
FEmatrices.SIS          = find((tab_region(:,3)+tab_region(:,2)) >= 1);%semi-infinite-space

labels_cell = get_labels([ext_plate_label,in_plate_label,BGPML_label],...
                         FILENAME);
FEmatrices.PlateExt  = find(labels_cell{1});
FEmatrices.PlateIn   = find(labels_cell{2});
FEmatrices.BGPML     = find(labels_cell{3});


tab_plate = [];
for ii=1:length(FEmatrices.plate_nodes)
    for jj=1:3
        tab_plate = [tab_plate 3*(FEmatrices.plate_nodes(ii)-1)+jj];
    end
end


plot3(FEmatrices.Nodes(FEmatrices.PlateExt,1),FEmatrices.Nodes(FEmatrices.PlateExt,2),FEmatrices.Nodes(FEmatrices.PlateExt,3),'+');
plot3(FEmatrices.Nodes(plate_nodes,1),FEmatrices.Nodes(plate_nodes,2),FEmatrices.Nodes(plate_nodes,3),'+');
plot3(FEmatrices.Nodes(FEmatrices.field,1),FEmatrices.Nodes(FEmatrices.field,2),FEmatrices.Nodes(FEmatrices.field,3),'+');

% indexing of the differents subspaces for partitionning
FEmatrices.indexu1        = 1:3:length(tab_plate);
FEmatrices.indexu2        = 2:3:length(tab_plate);
FEmatrices.indexu3        = 3:3:length(tab_plate);
FEmatrices.indexSIS       = (length(tab_plate)+1) : 1 : (length(tab_plate)+length(FEmatrices.SIS));
FEmatrices.indexBGPML     = zeros(length(FEmatrices.BGPML),1);
FEmatrices.indexBGfield    = zeros(length(FEmatrices.BG_nodes),1);

for ii=1:length(FEmatrices.indexBGfield)
    FEmatrices.indexBGfield(ii) = length(tab_plate)+find(FEmatrices.SIS==FEmatrices.BG_nodes(ii)); 
end

for ii=1:length(FEmatrices.BGPML)
    FEmatrices.indexBGPML(ii) = length(tab_plate)+find(FEmatrices.SIS==FEmatrices.BGPML(ii)); 
end


K = Kr(tab_plate,tab_plate) + 1i*Ki(tab_plate,tab_plate);
M = M(tab_plate,tab_plate);
FEmatrices.Hbg = Hbg(FEmatrices.BG_nodes,FEmatrices.BG_nodes);%needed for the calculation of the RHS
FEmatrices.Qbg = Qbg(FEmatrices.BG_nodes,FEmatrices.BG_nodes);%needed for the calculation of the RHS
Hpmlr = Hpmlr(FEmatrices.SIS,FEmatrices.SIS);
Hpmli = Hpmli(FEmatrices.SIS,FEmatrices.SIS);
Hpml = Hpmlr + 1i*Hpmli;
Qpmlr = Qpmlr(FEmatrices.SIS,FEmatrices.SIS);
Qpmli = Qpmli(FEmatrices.SIS,FEmatrices.SIS);
Qpml = Qpmlr + 1i*Qpmli;
C = sparse(size(M,1),size(Hpml,2));
C(FEmatrices.indexu3,:) = Ctmp(FEmatrices.plate_nodes,FEmatrices.SIS);
FEmatrices.Surf_matrix = Ctmp(FEmatrices.PlateIn, FEmatrices.PlateIn);


Kglob = sparse([K -C;...
                sparse(size(Hpml,1),size(K,2)) Hpml]);
Mglob = sparse([M sparse(size(M,1),size(Qpml,2));...
                param.rho*C' Qpml/param.c0^2]);

FEmatrices.LHS = {Kglob,Mglob};

FEmatrices.size_system = size(Kglob,1);

% FEmatrices.indexfield is the index of the nodes (of the RHS vector)  where is applied the BG
% pressure field
FEmatrices.field = FEmatrices.BG_nodes;
FEmatrices.indexfield = FEmatrices.indexBGfield;
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





