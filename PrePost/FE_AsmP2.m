function RHS = FE_AsmP2(FEmatrices,param,freq,theta,derivative_order)

[element_data, ~, Plan2D, normal_direction] = get_element_surface(FEmatrices, param);

gauss_point = [1/3. 1/3.];
gauss_coeff =  1./2.;

RHS_tmp = zeros(FEmatrices.size_system,1);
RHS = zeros(FEmatrices.size_system,1);

Area = 0;


for n=1:size(element_data,1)
    element_coordinates = [FEmatrices.Nodes(element_data(n,:), Plan2D(1)),...
                           FEmatrices.Nodes(element_data(n,:), Plan2D(2))];              
    [element_coordinates,idx] = sort_vertices(element_coordinates);
    element_data(n,:) = element_data(n,idx);
    detJ = 2*abs(polyarea(element_coordinates(:,1),...
                          element_coordinates(:,2)));
    Area = Area + detJ;
    
    for k=1:size(gauss_point,1)
        N = create_N(gauss_point(k,:), element_coordinates);
        X = N'*element_coordinates;
        for m=1:size(element_data,2) % 3 for P1 6 for P2
            RHS_tmp(element_data(n,m)) = RHS_tmp(element_data(n,m)) + ...
                                    FEmatrices.RHScoeffderiv_fun{1,derivative_order(1),...
                                                                   derivative_order(2)}(freq,...
                                                                                        theta,...
                                                                                        X(1),...
                                                                                        X(2))*...
                                    N(m)*...
                                    detJ*...
                                    gauss_coeff;
        end
    end
end

RHS(FEmatrices.indexPlateIn) = RHS_tmp(find(RHS_tmp));
       
end

function [element_data, element_center, plan2D, normal_direction] = get_element_surface(FEmatrices, param)
% list_nodes_surface: list of the nodes belonging to the surface
% element_data: 2D array (as much rows as elements, column1 = node1
%                                                   column2 = node2
%                                                   column3 = node3
connectivity_table = FEmatrices.connectivity+1;%FreeFem++ indexes nodes with 0 as first value (0+1=first value in Matlab)
list_nodes_surface = FEmatrices.PlateIn;
element_data = [];
table_elements_test = zeros(size(connectivity_table));
%connectivity_table = sort(connectivity_table,2);
for ii=1:size(connectivity_table,2)
    for jj=1:length(list_nodes_surface)
        idx = find(connectivity_table(:,ii)==list_nodes_surface(jj));
        table_elements_test(idx,ii) = list_nodes_surface(jj);
    end
end


for ii=1:size(table_elements_test,1)
    idx = find(table_elements_test(ii,:));
    if length(idx) == 6 % 6 nodes on the surface(for P2 triangles)
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
  
                           
normal_direction = find(indicator_Plan == min(indicator_Plan));

% %display part
% for ii=1:size(element_data,1)
%     X = FEmatrices.Nodes(element_data(ii,:),plan2D(1));
%     Y = FEmatrices.Nodes(element_data(ii,:),plan2D(2));
%     plot([X;element_center(ii,plan2D(1))], [Y;element_center(ii,plan2D(2))],'o');                         
%     hold on
% end    



end


function Ni = create_N(X,element_coordinates)

Ni = [(1-X(1)-X(2))*(1-2*X(1)-2*X(2));...
      4*X(1)*(1-X(1)-X(2));...
      X(1)*(2*X(1)-1);...
      4*X(1)*X(2);...
      X(2)*(2*X(2)-1);...
      4*X(2)*(1-X(1)-X(2))];
   
end


function [points, idx] = sort_vertices(points)

center = mean(points,1);
X = (points(:,1)-center(1));
Y = (points(:,2)-center(2));
sign_y = ones(size(points,1),1);
sign_y(find(Y<0)) = -1;
thetas = sign_y.*acos(X./(sqrt(X.^2+Y.^2)));
[thetas,idx] = sort(thetas);
points = points(idx,:);
[points, idx] = permute(points, idx);


% plot([-0.2,0.2],[0.2,-0.2],'+');
% hold on
% for ii=1:length(points)
%     plot(points(ii,1),points(ii,2),'o');
%     hold on
% end
end


function [P, idx] = permute(P, idx)

v1 = [P(2,1)-P(1,1);...
      P(2,2)-P(1,2)];
v2 = [P(3,1)-P(1,1);...
      P(3,2)-P(1,2)];

res = (v1(1)*v2(2) - v1(2)*v2(1))/norm(v1)/norm(v2);
tol = 1e-3;
if abs(res) > tol
    idx = idx([6 1 2 3 4 5]);
    P = P([6 1 2 3 4 5],:);
end

end



