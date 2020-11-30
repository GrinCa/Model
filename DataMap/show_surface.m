function show_surface(FEmatrices,param,SOLFE,nodes,Plan2D,normal_direction)

x = FEmatrices.Nodes(nodes,Plan2D(1));
y = FEmatrices.Nodes(nodes,Plan2D(2));

U_ortho = SOLFE(normal_direction:3:3*length(FEmatrices.plate_nodes),:,:);

f=figure;
for ii=1:param.nfreq
    plot3(x,y,real(U_ortho(nodes,ii,1)),"+");
    clf(f);
end



end