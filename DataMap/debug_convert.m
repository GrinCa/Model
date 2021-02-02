Nodes = load('Nodes_PML_Cube.txt');
mesh.file = 'Cube';
Pout = struct2cell(load('PoutWCAWE_Cube.mat'));
Pout = Pout{1};
convertGEO2VTK(mesh,Nodes,Pout,1);