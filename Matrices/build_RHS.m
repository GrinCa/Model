function RHS = build_RHS(f,theta,FEmatrices,RHScoeffderiv,derivative_order,direction)

RHS_function = RHScoeffderiv{1,derivative_order(1),derivative_order(2)};

Coordinates = FEmatrices.Nodes(FEmatrices.field,find(direction));
x1 = Coordinates(:,1);
x2 = Coordinates(:,2);


f = f*ones(size(FEmatrices.field,1),1);
theta = theta*ones(size(FEmatrices.field,1),1);

RHS = zeros(FEmatrices.size_system,1);
RHS(FEmatrices.indexfield) = RHS_function(f,theta,x1,x2);

end