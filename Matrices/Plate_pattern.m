function FEmatrices = Plate_pattern(FEmatrices,listLHS,param,FILENAME)

ndof = size(FEmatrices.Nodes,1);
% regions
acoustic_region = 1;
elastic_region = 2;
PML_region = 6;
%labels
embedding = 3;
coupling = 4;
extpressure = 5;


K     = listLHS{1}; % Stiffness matrix elastic domain
M     = listLHS{2}; % Mass matrix elastic domain
Hr    = listLHS{3};
%Hi    = listLHS{4};
Qr    = listLHS{4};
%Qi    = listLHS{6};
C1tmp = listLHS{5}; % coupling matrix
C2tmp = listLHS{6};


% get the arrays of the nodes belonging to a given region/label
tab_region = get_regions([acoustic_region,elastic_region,PML_region],ndof,FILENAME);
FEmatrices.acoustic_nodes     = find(tab_region(:,1)+tab_region(:,3));
FEmatrices.elastic_nodes      = find(tab_region(:,2));
FEmatrices.PML_nodes          = find(tab_region(:,3));

labels_cell = get_labels([embedding,coupling,extpressure],...
                          FILENAME);

FEmatrices.coupling_nodes     = find(labels_cell{2});
FEmatrices.extpressure_nodes  = find(labels_cell{3});

tab_plate = [];
for ii=1:length(FEmatrices.elastic_nodes)
    for jj=1:3
        tab_plate = [tab_plate 3*(FEmatrices.elastic_nodes(ii)-1)+jj];
    end
end

%for RHS calculation
FEmatrices.field = FEmatrices.extpressure_nodes;
FEmatrices.indexfield = zeros(length(FEmatrices.extpressure_nodes),1);
for ii=1:length(FEmatrices.extpressure_nodes)
    index_ii = find(FEmatrices.elastic_nodes == FEmatrices.extpressure_nodes(ii));
    FEmatrices.indexfield(ii) = 3*(index_ii-1)+1;
end


% FEmatrices.indexcoupling = zeros(length(FEmatrices.coupling_nodes),1);
% for ii=1:length(FEmatrices.coupling_nodes)
%     index_ii = find(FEmatrices.elastic_nodes == FEmatrices.coupling_nodes(ii));
%     FEmatrices.indexcoupling(ii) = 3*(index_ii-1)+1;
% end

%plot3(FEmatrices.Nodes(embedding_nodes,1),FEmatrices.Nodes(embedding_nodes,2),FEmatrices.Nodes(embedding_nodes,3),'+');

% indexing of the differents subspaces for partitionning
FEmatrices.indexu1   = 1:3:length(tab_plate);
FEmatrices.indexu2   = 2:3:length(tab_plate);
FEmatrices.indexu3   = 3:3:length(tab_plate);
FEmatrices.indexP    = (length(tab_plate)+1) : 1 : (length(tab_plate)+length(FEmatrices.acoustic_nodes));


K = K(tab_plate,tab_plate);
M = M(tab_plate,tab_plate);
Hr = Hr(FEmatrices.acoustic_nodes,FEmatrices.acoustic_nodes);
%Hi = Hi(FEmatrices.acoustic_nodes,FEmatrices.acoustic_nodes);
H = Hr;% + 1i*Hi;
Qr = Qr(FEmatrices.acoustic_nodes,FEmatrices.acoustic_nodes);
%Qi = Qi(FEmatrices.acoustic_nodes,FEmatrices.acoustic_nodes);
Q = Qr;% + 1i*Qi;
C = sparse(size(K,1),size(H,2));
C(FEmatrices.indexu1,:) = C2tmp(FEmatrices.elastic_nodes,FEmatrices.acoustic_nodes);

FEmatrices.SURFINT_1 = C1tmp(FEmatrices.extpressure_nodes, FEmatrices.extpressure_nodes);
FEmatrices.SURFINT_2 = C2tmp(FEmatrices.coupling_nodes, FEmatrices.coupling_nodes);

Kglob = sparse([K -C;...
                sparse(size(H,1),size(K,2)) H]);
Mglob = sparse([M sparse(size(M,1),size(Q,2));...
                param.rho*C' Q/param.c0^2]);

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



function FEmatrices = build_BGField(FEmatrices,param)
ndof = size(FEmatrices.Nodes,1);

% Cell of RHS
FEmatrices.RHS_BG = zeros(FEmatrices.size_system,param.nfreq,param.ntheta);

for ii=1:param.nfreq
    for jj=1:param.ntheta
        FEmatrices.RHS_BG(:,ii,jj) = build_RHS(param.freq(ii),param.theta(jj),FEmatrices,[1,1],param);
    end
end

end
