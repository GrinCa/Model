function RHS = fe_asm(FEmatrices,param,freq,theta,derivative_order)

[element_data, element_center, Plan2D] = get_element_surface(FEmatrices, param);

RHS = zeros(FEmatrices.size_system,1);

gauss_point = [1.0/3.0;...
               1.0/3.0];
gauss_coeff =  1.0/2.0;


for n=1:length(element_data)
    element_coordinates = [FEmatrices.Nodes(element_data(n,:), Plan2D(1)),...
                           FEmatrices.Nodes(element_data(n,:), Plan2D(2))];
    [N,detJ,element_coordinates] = create_N(gauss_point, element_coordinates);
    if detJ == 0
        disp("[WARNING] |J| <= 0")
    end
    X = N'*element_coordinates;
    for m=1:size(element_data,2)
        RHS(element_data(n,m)) = RHS(element_data(n,m)) + ...
                                 FEmatrices.RHScoeffderiv_fun{1,derivative_order(1),derivative_order(2)}(freq,...
                                                                                                         theta,...
                                                                                                         X(1),...
                                                                                                         X(2))*...
                                 N(m)*...
                                 detJ*...
                                 gauss_coeff;
    end
end
           
end

function [element_data, element_center, plan2D] = get_element_surface(FEmatrices, param)
% list_nodes_surface: list of the nodes belonging to the surface
% element_data: 2D array (as much rows as elements, column1 = node1
%                                                   column2 = node2
%                                                   column3 = node3
connectivity_table = FEmatrices.connectivity+1;%FreeFem++ indexes nodes with 0 as first value (0+1=first value in Matlab)
list_nodes_surface = FEmatrices.PlateIn;
element_data = [];
table_elements_test = zeros(size(connectivity_table));
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


% %Calculation of the surface of the 2D triangle of the surface
indicator_Plan = [norm(FEmatrices.Nodes(FEmatrices.PlateIn,1)-FEmatrices.Nodes(FEmatrices.PlateIn(1),1)),...
                  norm(FEmatrices.Nodes(FEmatrices.PlateIn,2)-FEmatrices.Nodes(FEmatrices.PlateIn(1),2)),...
                  norm(FEmatrices.Nodes(FEmatrices.PlateIn,3)-FEmatrices.Nodes(FEmatrices.PlateIn(1),3))];
plan2D = find(indicator_Plan>min(indicator_Plan)); % determines whether the plan is (x,y) (x,z) (y,z) ...
                           % |_>> min(indicator_Plan) = 0 if there is no "numerical epsilon"
  
                           
%normal_direction = find(indicator_Plan == min(indicator_Plan));

% %display part
% for ii=1:size(element_data,1)
%     X = FEmatrices.Nodes(element_data(ii,:),plan2D(1));
%     Y = FEmatrices.Nodes(element_data(ii,:),plan2D(2));
%     plot([X;element_center(ii,plan2D(1))], [Y;element_center(ii,plan2D(2))],'o');                         
%     hold on
% end    

end


function [Ni,detJ,element_coordinates] = create_N(X,element_coordinates)


% L = [1                        1                        1                       ;...
%      element_coordinates(1,1) element_coordinates(2,1) element_coordinates(3,1);...
%      element_coordinates(1,2) element_coordinates(2,2) element_coordinates(3,2)];
%  
% uno = ones(1,size(X,2));
% X = [uno;X];
% 
% Ni = L\X;

Ni = [1-X(1)-X(2);...
      X(1);...
      X(2)];

detJ = (element_coordinates(2,1)-element_coordinates(1,1))*(element_coordinates(3,2)-element_coordinates(1,2))-...
       (element_coordinates(3,1)-element_coordinates(1,1))*(element_coordinates(2,2)-element_coordinates(1,2));

if detJ <= 0
    element_coordinates = element_coordinates([1 3 2],:);
    detJ = (element_coordinates(2,1)-element_coordinates(1,1))*(element_coordinates(3,2)-element_coordinates(1,2))-...
           (element_coordinates(3,1)-element_coordinates(1,1))*(element_coordinates(2,2)-element_coordinates(1,2));
end
   
end
