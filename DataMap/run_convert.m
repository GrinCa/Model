function run_convert(param, mesh, flag)
    clear FEmatrices;
    sizemesh_ARRAY = load('sizemesh.txt');
    sizemesh = sizemesh_ARRAY(end);
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/','DATA_sizemesh_',num2str(sizemesh),'.mat']));
    FEmatrices = DATA{1};
    param = DATA{2};
    vtk_config = str2func([mesh.file '_vtk']);
    param = vtk_config(FEmatrices,param);
    if flag.calculateFE
        SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLFE = SOLFE{1};      
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLFE,param.VTK.PARTITION,param,param.VTK.range);
    end
end