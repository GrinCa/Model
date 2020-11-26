function [RHS,BG_field] = build_RHS(f,theta,FEmatrices,derivative_order,param)

RHS_function = FEmatrices.RHScoeffderiv_fun{1,derivative_order(1),derivative_order(2)};


Coordinates = FEmatrices.Nodes(FEmatrices.field,find(param.direction));


x1 = Coordinates(:,1);
x2 = Coordinates(:,2);

f_array = f*ones(length(FEmatrices.field),1);
theta_array = theta*ones(length(FEmatrices.field),1);

RHS = zeros(FEmatrices.size_system,1);
Z = FEmatrices.Hbg - (2*pi*f/param.c0)^2*FEmatrices.Qbg;
RHStmp = RHS_function(f_array,theta_array,x1,x2);
%RHS(FEmatrices.indexfield) = 2*RHStmp.*FEmatrices.RHS;
RHS(FEmatrices.indexfield) = -Z*RHS_function(f_array,theta_array,x1,x2);
RHS(FEmatrices.indexBGPML) = 0;

BG_field = RHS_function(f_array,theta_array,x1,x2);
end