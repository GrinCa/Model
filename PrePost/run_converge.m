function run_converge(param)

    clear FEmatrices SOLFE SOLMDWCAWE;
    sizemesh_ARRAY = load('sizemesh.txt');
    fid1 = fopen('converge1.txt','wt');
    fid2 = fopen('converge2.txt','wt');
    if flag.calculateFE == 1
        for ii=1:length(sizemesh_ARRAY)
            sizemesh = sizemesh_ARRAY(ii);    
            SOL = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLFE_',mesh.file,'_sizemesh_',num2str(sizemesh),'.mat']));
            SOL = SOL{1};
        end
    elseif flag.calculateMDWCAWE
        for ii=1:length(sizemesh_ARRAY)
            sizemesh = sizemesh_ARRAY(ii);
            SOL = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLMDWCAWE_',mesh.file,'_nvec_',num2str(param.n_sub_range*nvecfreq),'_sizemesh_',num2str(sizemesh),'.mat']));
            SOL = SOL{1};
        end
    end
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/','DATA_sizemesh_',num2str(sizemesh),'.mat']));
    FEmatrices = DATA{1};
    param = DATA{2};
    ndof = size(FEmatrices.Nodes,1);
    ndof_acoustic = length(FEmatrices.indexp);
    %ndof_acoustic = length(FEmatrices.acoustic_nodes);
    acoustic_volume = load('acoustic_volume.txt');
    Pressure = zeros(ndof,param.nfreq);
    Pressure(FEmatrices.acoustic_nodes,:) = real(SOL(FEmatrices.indexp,:));
    MQP1 = Pressure(FEmatrices.acoustic_nodes,:)'*FEmatrices.Q*Pressure(FEmatrices.acoustic_nodes,:)/((4e-10)*acoustic_volume);
    MQP1 = 10*log10(MQP1);
    fprintf(fid1,strcat(num2str(ndof),'\t',num2str(real(diag(MQP1)')),'\n'));
    MQP2 = Pressure'*Pressure/4/(10e-10)/2/ndof_acoustic;
    MQP2 = 10*log10(MQP2);
    fprintf(fid2,strcat(num2str(ndof),'\t',num2str(real(diag(MQP2)')),'\n'));
    fclose(fid1);
    fclose(fid2);

end

%     sizemesh_file = load('sizemesh.txt');
%     meanFE = cell(length(sizemesh_file),1);
%     title_VALUES_1 = cell(length(sizemesh_file),1);
%     argcomp.label = cell(length(sizemesh_file),1);
%     for ii=1:length(sizemesh_file)
%         DATA = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/DATA_sizemesh_',num2str(sizemesh_file(ii)),'.mat']));
%         FEmatrices = DATA{1};
%         param = DATA{2};
%         SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/SOLFE_sizemesh_',num2str(sizemesh_file(ii)),'.mat']));
%         SOLFE = SOLFE{1};
%         meanFE{ii} = mean(real(SOLFE(FEmatrices.PlateCavity_nodes,:)),1);
%         argcomp.label{ii} = [num2str(size(FEmatrices.Nodes,1)) ' ndofs'];
%     end
% 
%     argcomp.type = 'converge';
%     argcomp.xlabel = 'freq';
%     argcomp.ylabel = 'mean pressure (Pa)';
%     argcomp.title = '';
%     argcomp.name_plot = ['Convergence_FE_[' replace(num2str(sizemesh_file'),' ','_') ']'];
%     argcomp.external_plot.is_needed = false;
%     show_graph(argcomp,meanFE,mesh,param);