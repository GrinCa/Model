function FEmatrices = Modelv3_pattern(FEmatrices,listLHS,param,FILENAME)

ndof = size(FEmatrices.Nodes,1);
% regions
plate_region = 1;
%labels
in_plate_label = 2;
ext_plate_label = 3;
embedding = 4;

K     = listLHS{1}; % Stiffness matrix elastic domain
M     = listLHS{2}; % Mass matrix elastic domain
C     = listLHS{3}; % Coupling matrix for surface calculation (see get_power.m)

%connectivity
FEmatrices.connectivity = load(['Matrices/',FILENAME,'/connectivity_table.txt']);


% get the arrays of the nodes belonging to a given region/label
tab_region = get_regions([plate_region],ndof,FILENAME);
FEmatrices.plate_nodes  = find(tab_region(:,1));

labels_cell = get_labels([in_plate_label,ext_plate_label,embedding],...
                         FILENAME);
FEmatrices.PlateIn   = find(labels_cell{1});
FEmatrices.PlateExt  = find(labels_cell{2});
FEmatrices.embedding = find(labels_cell{3});


tab_plate = [];
for ii=1:length(FEmatrices.plate_nodes)
    for jj=1:3
        tab_plate = [tab_plate 3*(FEmatrices.plate_nodes(ii)-1)+jj];
    end
end


% plot3(FEmatrices.Nodes(FEmatrices.PlateExt,1),FEmatrices.Nodes(FEmatrices.PlateExt,2),FEmatrices.Nodes(FEmatrices.PlateExt,3),'+');
% plot3(FEmatrices.Nodes(plate_nodes,1),FEmatrices.Nodes(plate_nodes,2),FEmatrices.Nodes(plate_nodes,3),'+');
% plot3(FEmatrices.Nodes(FEmatrices.field,1),FEmatrices.Nodes(FEmatrices.field,2),FEmatrices.Nodes(FEmatrices.field,3),'+');

% indexing of the differents subspaces for partitionning
FEmatrices.indexu1        = 1:3:length(tab_plate);
FEmatrices.indexu2        = 2:3:length(tab_plate);
FEmatrices.indexu3        = 3:3:length(tab_plate);

FEmatrices.Surf_matrix = C(FEmatrices.PlateIn, FEmatrices.PlateIn);


Kglob = K(tab_plate,tab_plate);
Mglob = M(tab_plate,tab_plate);

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





