function [rayleigh_matrices, element_data] = RIM(FEmatrices,param)
% rayleigh_matrices: cell of matrices (as many as frequencies)
rayleigh_matrices = cell(param.nfreq,1);
% normal_direction: integer (1,2 or 3) whether the normal of the plate is 
% x,y or z

% get all data needed for the calculation of the matrix
[element_data, element_center] = get_element_surface(FEmatrices, param);

% for ii=1:size(element_data)
%     plot3(FEmatrices.Nodes(element_data(ii,:),1),FEmatrices.Nodes(element_data(ii,:),2),FEmatrices.Nodes(element_data(ii,:),3),'o');
%     hold on
% end

n_elem = length(element_data);

unit = ones(size(FEmatrices.SurfExt_matrix,1),1);
Se = abs( unit'*FEmatrices.SurfExt_matrix*unit / n_elem );

disp("Start Rayleigh calculation");

for n=1:param.nfreq
    rayleigh_matrices{n} = zeros(n_elem, n_elem);
    omega = 2*pi*param.freq(n);
    k = omega/param.c0;
    coeff = Se^2*omega^2*param.rho0 / (4*pi*param.c0);
    disp(['Rayleigh Matrix [f = ',num2str(param.freq(n)),' Hz]']);
    for ii=1:n_elem
        for jj=1:(ii-1)
            rij = norm(element_center(ii,:) - element_center(jj,:));
            rayleigh_matrices{n}(ii,jj) = sin(k*rij)/k/rij;
        end
    end
    rayleigh_matrices{n} = coeff * ( rayleigh_matrices{n} + ...
                                     rayleigh_matrices{n}'+ ...
                                     eye(size(rayleigh_matrices{n},1)) );
end

end

function [element_data, element_center] = get_element_surface(FEmatrices, param)
% list_nodes_surface: list of the nodes belonging to the surface
% element_data: 2D array (as much rows as elements, column1 = number of the
%                                                             element
%                                                   column2 = node1
%                                                   column3 = node2
%                                                   column4 = node3
connectivity_table = FEmatrices.connectivity+1;%FreeFem++ indexes nodes with 0 as first value (0+1=first value in Matlab)
list_nodes_surface = FEmatrices.PlateExt;
element_data = [];
table_elements_test = zeros(size(connectivity_table,1),...
                            size(connectivity_table,2));
connectivity_table = sort(connectivity_table,2);
for ii=1:length(list_nodes_surface)
    for jj=1:size(connectivity_table,2)
        idx = find(connectivity_table(:,jj)==list_nodes_surface(ii));
        table_elements_test(idx,jj) = list_nodes_surface(ii);
    end
end


for ii=1:size(table_elements_test,1)
    idx = find(table_elements_test(ii,:));
    if length(idx) == 3 % 3 nodes on the surface(for P1 triangles)
        element_data = [element_data;table_elements_test(ii,idx)];
    end
end



% Calculation of the center of each element
element_center = zeros(size(element_data,1),3);
for ii=1:size(element_data,1)
    element_center(ii,:) = mean(FEmatrices.Nodes(element_data(ii,:),:),1);
end

% 
% % %Calculation of the surface of the 2D triangle of the surface
% indicator_Plan = [norm(FEmatrices.Nodes(FEmatrices.PlateExt,1)-FEmatrices.Nodes(FEmatrices.PlateExt(1),1)),...
%                   norm(FEmatrices.Nodes(FEmatrices.PlateExt,2)-FEmatrices.Nodes(FEmatrices.PlateExt(1),2)),...
%                   norm(FEmatrices.Nodes(FEmatrices.PlateExt,3)-FEmatrices.Nodes(FEmatrices.PlateExt(1),3))];
% Plan2D = find(indicator_Plan>min(indicator_Plan)); % determines whether the plan is (x,y) (x,z) (y,z) ...
%                            % |_>> min(indicator_Plan) = 0 if there is no "numerical epsilon"
%  


%display part
% for ii=1:size(element_data,1)
%     X = FEmatrices.Nodes(element_data(ii,:),plan2D(1));
%     Y = FEmatrices.Nodes(element_data(ii,:),plan2D(2));
%     plot([X;element_center(ii,plan2D(1))], [Y;element_center(ii,plan2D(2))],'o');                         
%     hold on
% end    

end









