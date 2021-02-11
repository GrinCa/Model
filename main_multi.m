%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Test case using WCAWE                             %
%                                                                         %
%                             March 2020                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% Init main program
%--------------------------------------------------------------------------


% clear FEmatrices param flag;
% warning('off', 'MATLAB:nearlySingularMatrix');

%--------------------------------------------------------------------------
% Folders
%--------------------------------------------------------------------------

% Add folders for Mumps, WCAWE, and Mesh functions/files

Meshfold = 'Matrices';
WCAWEfold = 'WCAWE';
Mumps = 'Mumps';
DataMap = 'DataMap';
Derivatives = 'Derivatives';
STL = 'STL';
Config = 'Config';

addpath(genpath(strcat(pwd,'/',Meshfold)));
addpath(genpath(strcat(pwd,'/',WCAWEfold)));
addpath(genpath(strcat(pwd,'/',Mumps)));
addpath(genpath(strcat(pwd,'/',DataMap)));
addpath(genpath(strcat(pwd,'/',Derivatives)));
addpath(genpath(strcat(pwd,'/',STL)));
addpath(genpath(strcat(pwd,'/',Config)));


%--------------------------------------------------------------------------
% Input data for the problem
%--------------------------------------------------------------------------
% Input files for mesh and FE matrices

sizemesh = load('sizemesh.txt');
sizemesh = sizemesh(end);
config = str2func([mesh.file '_config']);
[flag, param] = config();


% generation of the different folder to store data if they don't already exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param = genfolders(mesh,param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Build intervals for the Parametric sweep
%--------------------------------------------------------------------------
param = build_interval(param);

%-----------------------------------------------------------------------------------------
[FEmatrices.LHScoeffderiv_fun, FEmatrices.RHScoeffderiv_fun] = get_coeff_deriv_matrices();
%-----------------------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Matrices calculated with Freefem++
%--------------------------------------------------------------------------
if flag.getmatrices
    [FEmatrices,ndof,flag] = get_matrices(FEmatrices,flag,mesh,param);
    Nodes = FEmatrices.Nodes;
    LHS = FEmatrices.LHS;
    nLHS = length(LHS);
    param.idx_out = 1:1:FEmatrices.size_system;
    % Save data only for FE solution
    save(['Matrices/',mesh.file,'/',param.path1,'/','DATA_sizemesh_',num2str(sizemesh),'.mat'],'FEmatrices','param');
end

if flag.recalculated

    %--------------------------------------------------------------------------
    % Calculate reference finite element parametric sweep
    %--------------------------------------------------------------------------
    
    [FEmatrices.LHScoeffderiv_fun, FEmatrices.RHScoeffderiv_fun] = get_coeff_deriv_matrices();
    BG_field = zeros(size(FEmatrices.Nodes,1),param.nfreq,param.ntheta);

    if flag.calculateFE == 1
       disp('**************************');
       disp('* Compute FE calculation *');
       disp('**************************');
       SOLFE = zeros(FEmatrices.size_system,param.nfreq,param.ntheta); %size ndof,nfreq
       % Parametric loop calculation
       id = initmumps;
       id = zmumps(id);
       id.JOB = 1; %analysis
       matrix_analysis = LHS{1} + LHS{2};
       id = zmumps(id,matrix_analysis); % perform analysis
       for ii=1:param.nfreq
           for jj=1:param.ntheta
               disp(['[FE] [Frequency, Theta] = [',num2str(param.freq(ii)),',',num2str(param.theta(jj)/pi*180),']']);
               Aglob = sparse(size(LHS{1},1),size(LHS{1},2));
               for kk = 1:nLHS
                  Aglob = Aglob + param.coeff_LHS{kk}(param.freq(ii))*LHS{kk};
               end
               id.JOB = 5;
               id.RHS = FE_AsmP2(FEmatrices,param,param.freq(ii),param.theta(jj),[1 1]);
               id = zmumps(id,Aglob);
               resp_P = id.SOL;
               SOLFE(:,ii,jj) = resp_P;
           end
       end
       id.JOB = -2; id = zmumps(id); %delete the instance of Mumps
    end


    %--------------------------------------------------------------------------
    % Calculate WCAWE parametric sweep
    %--------------------------------------------------------------------------

    if flag.calculateMDWCAWE || flag.calculateWCAWE

        % Initialize array of RHSderiv, Cell array linear comb coefficients derivative functions
        deriv_deg = [param.nvecfreqmax, param.nvecthetamax];

        if exist('Derivatives/derivative_orders.mat','file') ~= 2
            disp('*************************************');
            disp('* Recalculate all cross derivatives *');
            disp('*************************************');
            create_cross_derivatives(nLHS,...
                                     param.coeff_LHS,...
                                     param.coeff_RHS,deriv_deg,...
                                     {'f';'theta'},{'f,theta';...
                                                    'f,theta';...
                                                    'f,theta,x1,x2'});
        end
        load('Derivatives/derivative_orders.mat');
        if ~isempty(find(derivative_orders-deriv_deg<0))
            disp('*************************************');
            disp('* Recalculate all cross derivatives *');
            disp('*************************************');
            create_cross_derivatives(LHS,...
                                     param.coeff_LHS,...
                                     param.coeff_RHS,deriv_deg,...
                                     {'f';'theta'},{'f,theta';...
                                                    'f,theta';...
                                                    'f,theta,x1,x2'});
        end

        for nvecfreq=param.vecfreqrange
            for nvectheta=param.vecthetarange
                %%%%%
                param.nvecfreq = nvecfreq;
                param.nvectheta = nvectheta;
                %%%%%
                [FEmatrices.LHScoeffderiv_fun, FEmatrices.RHScoeffderiv_fun] = get_coeff_deriv_matrices();
                [FEmatrices.LHScoeffderiv, FEmatrices.RHSderiv] = fill_array_WCAWE(FEmatrices,param);

                %-----------------------------------------------------------------------
                %MDWCAWE Basis Julien
                %-----------------------------------------------------------------------
                if flag.calculateMDWCAWE
                    disp('****************************************');
                    disp('* Compute calculation of MDWCAWE basis *');
                    disp('****************************************');
                    arg.algo = 'MDWCAWE';
                    SOLMDWCAWE = zeros(length(param.idx_out),length(param.freq),length(param.theta));

                    for n=1:length(param.interval_construct{1}) %loop over frequency intervals
                        for m=1:length(param.interval_construct{2}) % loop over angle intervals
                            MDWCAWEtmp = build_basis(FEmatrices,param,...
                                                     nvecfreq,nvectheta,...
                                                     param.interval_construct{1}{n},param.interval_construct{2}{m},...
                                                     arg); 
                            disp('**** Compute MDWCAWE projection ****');
                            SOLMDWCAWE(:,param.interval_index{1}{n},param.interval_index{2}{m}) = Solve_MDWCAWE(FEmatrices,param,MDWCAWEtmp,param.sub_interval{1}{n},param.sub_interval{2}{m});
                            disp('**** MDWCAWE projection done ****');
                        end
                    end
                end

                %-----------------------------------------------------------------------
                %WCAWE Basis Romain
                %-----------------------------------------------------------------------
                if flag.calculateWCAWE
                    disp('**************************************');
                    disp('* Compute calculation of WCAWE basis *');
                    disp('**************************************');
                    % Setup list of derivative coefficients and RHS derivatives
                    arg.algo = 'WCAWE';
                    SOLWCAWE = zeros(length(param.idx_out),length(param.freq),length(param.theta));

                    for n=1:length(param.interval_construct{1}) %loop over frequency intervals
                        for m=1:length(param.interval_construct{2}) % loop over angle intervals
                            WCAWEtmp = build_basis(FEmatrices,param,...
                                                   nvecfreq,nvectheta,...
                                                   param.interval_construct{1}{n},param.interval_construct{2}{m},...
                                                   arg); 
                            disp('**** Compute MDWCAWE projection ****');
                            SOLWCAWE(:,param.interval_index{1}{n},param.interval_index{2}{m}) = Solve_MDWCAWE(FEmatrices,param,WCAWEtmp,param.sub_interval{1}{n},param.sub_interval{2}{m});
                            disp('**** MDWCAWE projection done ****');
                        end
                    end
                end
            end
       end
    end

    %--------------------------------------------------------------------------
    % Saves
    %--------------------------------------------------------------------------
    % FE solution
    if flag.calculateFE
        save(['Matrices/',mesh.file,'/',param.path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat'],'SOLFE');
    end

    if flag.calculateMDWCAWE
        id_sample = ['SOLMDWCAWE_sizemesh_',num2str(sizemesh),'_',param.interval_detail_str '.mat'];
        for ii=1:length(param.vecfreqrange)
            for jj=1:length(param.vecthetarange)
                save(['Matrices/',mesh.file,'/',param.path1,'/[',replace(num2str(param.freqref),' ','_'),'][',replace(num2str(int16(180/pi*param.thetaref)),' ','_'),']/[',num2str(param.vecfreqrange(ii)),'_',num2str(param.vecthetarange(jj)),']/',id_sample], 'SOLMDWCAWE');
            end
        end
    end

    if flag.calculateWCAWE
        id_sample = ['SOLWCAWE_sizemesh_',num2str(sizemesh),'_',param.interval_detail_str '.mat'];
        for ii=1:length(param.vecfreqrange)
            for jj=1:length(param.vecthetarange)
                save(['Matrices/',mesh.file,'/',param.path1,'/[',replace(num2str(param.freqref),' ','_'),'][',replace(num2str(int16(180/pi*param.thetaref)),' ','_'),']/[',num2str(param.vecfreqrange(ii)),'_',num2str(param.vecthetarange(jj)),']/',id_sample], 'SOLWCAWE');
            end
        end
    end

    % save FEmatrices which contains all the data of the simulation for
    % each mesh

end

%--------------------------------------------------------------------------
% Post processing
%--------------------------------------------------------------------------

%show_surface(FEmatrices,param,SOLFE,FEmatrices.PlateExt,[1 2],3);

if flag.converge
    clear FEmatrices SOLFE SOLMDWCAWE;
    sizemesh_ARRAY = load('sizemesh.txt');
    fid1 = fopen('converge1.txt','wt');
    fid2 = fopen('converge2.txt','wt');
    if flag.calculateFE == 1
        for ii=1:length(sizemesh_ARRAY)
            sizemeshtmp = sizemesh_ARRAY(ii);    
            SOL = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLFE_',mesh.file,'_sizemesh_',num2str(sizemesh),'.mat']));
            SOL = SOL{1};
        end
    elseif flag.calculateMDWCAWE
        for ii=1:length(sizemesh_ARRAY)
            sizemeshtmp = sizemesh_ARRAY(ii);
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

if flag.plotMQP
    fid = fopen('converge.txt','rt');
    while true
        line = fgets(fid);
        if line == -1
            break;
        else
            line = str2num(strtrim(line));
            plot(param.freq,line(2:end),'DisplayName',['ndof = ' num2str(line(1))]);
            hold on
        end
    end
    
    xlabel('Frequency (Hz)');
    ylabel('Mean quadratic pressure (dB)');
    legend();
    hold off
    fclose(fid);
end


%--------------------------------------------------------------------------
% convert
%--------------------------------------------------------------------------
if flag.convert2VTK
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

% TL = get_STL(FEmatrices,param,SOLFE);
% plot(param.freq,TL);

if flag.calculateTL
    clear FEmatrices SOLFE SOLWCAWE SOLMDWCAWE;
    arg.sizemesh = sizemesh;
    arg.type = 'preload';
    FEmatrices = IO_data(arg,param,mesh);
    
    VALUES_SOL = cell(2,2); % index 2
    VALUES_name = cell(1,2);
    % is for MDWCAWE and WCAWE
    
    try
        arg.type = 'load_FE';
        SOLFE = IO_data(arg,param,mesh);
        arg.type = 'calculateTL';
        TL_FE = post_process(FEmatrices,param,arg,SOLFE);
        VALUES_SOL{1,1} = TL_FE;
        VALUES_SOL{2,1} = TL_FE;
        VALUES_name{1} = ['FE'];
    catch
        disp("SOLFE not find for the given parameters");
        TL_FE = zeros(param.nfreq,param.ntheta);
        VALUES_SOL{1,1} = TL_FE;
        VALUES_SOL{2,1} = TL_FE;
        VALUES_name{1} = ['FE'];
    end
    counter = 2;
    try
        arg.type = 'load_MDWCAWE';
        SOLMDWCAWE = IO_data(arg,param,mesh);
        arg.type = 'calculateTL';
        for ii=1:length(param.vecfreqrange)
            for jj=1:length(param.vecthetarange)
                TL_MDWCAWE = post_process(FEmatrices,param,arg,SOLMDWCAWE{ii,jj});
                VALUES_SOL{1,counter} = TL_MDWCAWE;
                VALUES_name{counter} = ['FPS [' num2str(param.vecfreqrange(ii)) num2str(param.vecthetarange(jj)) ']'];
                counter = counter +1;
            end
        end
    catch
        disp("SOLMDWCAWE not find for the given parameters");
        for ii=1:length(param.vecfreqrange)
            for jj=1:length(param.vecthetarange)
            TL_MDWCAWE = zeros(param.nfreq,param.ntheta);
            VALUES_SOL{1,counter} = TL_MDWCAWE;
            VALUES_name{counter} = ['FPS [' num2str(param.vecfreqrange(ii)) num2str(param.vecthetarange(jj)) ']'];
            counter = counter +1;
            end
        end
    end
    counter = 2;
    try
        for ii=1:length(param.vecfreqrange)
            for jj=1:length(param.vecthetarange)
                arg.type = 'load_WCAWE';
                SOLWCAWE = IO_data(arg,param,mesh);
                arg.type = 'calculateTL';
                VALUES_SOL{2,counter} = post_process(FEmatrices,param,arg,SOLWCAWE{ii,jj});
                counter = counter +1;
            end
        end
    catch
        disp("SOLWCAWE not find for the given parameters");
        for ii=1:length(param.vecfreqrange)
            for jj=1:length(param.vecthetarange)
                TL_WCAWE = zeros(param.nfreq,param.ntheta);
                VALUES_SOL{2,counter} = TL_WCAWE;
                counter = counter +1;
            end
        end
    end
   
    argcomp.ylabel = 'TL';
    argcomp.title = 'Comparison FE MDWCAWE';
    argcomp.type = 'plotTL';
    argcomp.split = 0;
    argcomp.name_plot = 'Comparison_FE_MDWCAWE';
    argcomp.label = VALUES_name;
    argcomp.external_plot.is_needed = false;
    show_graph(argcomp,VALUES_SOL(1,:),mesh,param);
    
    argcomp.ylabel = 'TL';
    argcomp.title = 'Comparison FE WCAWE';
    argcomp.type = 'plotTL';
    argcomp.split = 0;
    argcomp.name_plot = 'Comparison_FE_WCAWE';
    argcomp.label = VALUES_name;
    argcomp.external_plot.is_needed = false;
    show_graph(argcomp,VALUES_SOL(2,:),mesh,param);
    
    argtmp.type = 'rel_error';
    argtmp.REF_SOLUTION = VALUES_SOL(1);
    argtmp.APPROX_SOLUTION = VALUES_SOL(1,2:end);
    argcomp.type = 'plotTL';
    argcomp.ylabel = 'relative error';
    argcomp.title = '';
    argcomp.split = 1;
    argcomp.name_plot = 'Relative_error';
    argcomp.label = VALUES_name(2:end);
    argcomp.external_plot.is_needed = false;
    show_graph(argcomp,post_process(FEmatrices,param,argtmp),mesh,param);
end


if flag.eigen
    eigs(LHS{1},LHS{2},10,'sm')
end


if flag.converge_sizemesh
    sizemesh_file = load('sizemesh.txt');
    meanFE = cell(length(sizemesh_file),1);
    title_VALUES_1 = cell(length(sizemesh_file),1);
    argcomp.label = cell(length(sizemesh_file),1);
    for ii=1:length(sizemesh_file)
        DATA = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/DATA_sizemesh_',num2str(sizemesh_file(ii)),'.mat']));
        FEmatrices = DATA{1};
        param = DATA{2};
        SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/SOLFE_sizemesh_',num2str(sizemesh_file(ii)),'.mat']));
        SOLFE = SOLFE{1};
        meanFE{ii} = mean(real(SOLFE(FEmatrices.PlateCavity_nodes,:)),1);
        argcomp.label{ii} = [num2str(size(FEmatrices.Nodes,1)) ' ndofs'];
    end

    argcomp.type = 'converge';
    argcomp.xlabel = 'freq';
    argcomp.ylabel = 'mean pressure (Pa)';
    argcomp.title = '';
    argcomp.name_plot = ['Convergence_FE_[' replace(num2str(sizemesh_file'),' ','_') ']'];
    argcomp.external_plot.is_needed = false;
    show_graph(argcomp,meanFE,mesh,param);
end

if flag.normalized_error
    VALUES = cell(1,2);
    arg.sizemesh = sizemesh;
    arg.type = 'preload';
    FEmatrices = IO_data(arg,param,mesh);
    arg.type = 'load_FE';
    arg.REF_SOLUTION = IO_data(arg,param,mesh);
    arg.type = 'load_MDWCAWE';
    arg.APPROX_SOLUTION = IO_data(arg,param,mesh);
    arg.type = 'normalize_error';
    VALUES{1} = log10(post_process(FEmatrices,param,arg));
    arg.split = false;
    arg.label = {'MDWCAWE','WCAWE'};
    arg.zlabel = 'error';
    arg.title = {'FE/MDWCAWE','FE/WCAWE'};
    arg.external_plot.is_needed = false;
    arg.name_plot = 'relative error norm';
    %show_graph(arg,VALUES,mesh,param);
    
    
    arg.type = 'load_WCAWE';
    arg.APPROX_SOLUTION = IO_data(arg,param,mesh);
    arg.type = 'normalize_error';
    arg.split = true;
    VALUES{2} = log10(post_process(FEmatrices,param,arg));
    show_graph(arg,VALUES,mesh,param);
end













