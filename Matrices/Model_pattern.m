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





