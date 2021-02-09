function param = Modelv4_vtk(FEmatrices, param)

% convert2vtk

param.VTK.PARTITION = cell(1);

param.VTK.PARTITION{1} = {'vector',...
                          [FEmatrices.indexu1;...
                           FEmatrices.indexu2;...
                           FEmatrices.indexu3],...
                           FEmatrices.plate_nodes,...
                          'DISPLACEMENT'};
 
param.VTK.range = {1:1:1, 1:1:1};



end