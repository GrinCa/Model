function data = IO_data(arg,param,mesh)

if strcmp(arg.type,'load_FE')
    data = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/SOLFE_sizemesh_',num2str(arg.sizemesh),'.mat']));
    data = data{1};
end

if strcmp(arg.type,'load_MDWCAWE')
    id_sample = ['SOLMDWCAWE_sizemesh_',num2str(arg.sizemesh),'_',param.interval_detail_str,'.mat'];
    data = cell(length(param.vecfreqrange),length(param.vecthetarange));
    for ii=1:length(param.vecfreqrange)
        for jj=1:length(param.vecthetarange)
            SOLtmp = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/[',replace(num2str(param.freqref),' ','_'),'][',replace(num2str(int16(180/pi*param.thetaref)),' ','_'),']/[',num2str(param.vecfreqrange(ii)),'_',num2str(param.vecthetarange(jj)),']/',id_sample]));
            data{ii,jj} = SOLtmp{1};
        end
    end
end

if strcmp(arg.type,'load_WCAWE')
    id_sample = ['SOLWCAWE_sizemesh_',num2str(arg.sizemesh),'_',param.interval_detail_str,'.mat'];
    data = cell(length(param.vecfreqrange),length(param.vecthetarange));
    for ii=1:length(param.vecfreqrange)
        for jj=1:length(param.vecthetarange)
            SOLtmp = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/[',replace(num2str(param.freqref),' ','_'),'][',replace(num2str(int16(180/pi*param.thetaref)),' ','_'),']/[',num2str(param.vecfreqrange(ii)),'_',num2str(param.vecthetarange(jj)),']/',id_sample]));
            data{ii,jj} = SOLtmp{1};
        end
    end
end

if strcmp(arg.type,'preload')
    data = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/','DATA_sizemesh_',num2str(arg.sizemesh),'.mat']));
end

if strcmp(arg.type,'save_plot')
    saveas(arg.gcf,['Matrices/',mesh.file,'/',param.path1,'/[',replace(num2str(param.freqref),' ','_'),'][',replace(num2str(int16(180/pi*param.thetaref)),' ','_'),']/',arg.save_name]);
end

end


