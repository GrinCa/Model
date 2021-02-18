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
Timing = 'Timing';
PrePost = 'PrePost';

addpath(genpath(strcat(pwd,'/',Meshfold)));
addpath(genpath(strcat(pwd,'/',WCAWEfold)));
addpath(genpath(strcat(pwd,'/',Mumps)));
addpath(genpath(strcat(pwd,'/',DataMap)));
addpath(genpath(strcat(pwd,'/',Derivatives)));
addpath(genpath(strcat(pwd,'/',STL)));
addpath(genpath(strcat(pwd,'/',Config)));
addpath(genpath(strcat(pwd,'/',Timing)));
addpath(genpath(strcat(pwd,'/',PrePost)));

%--------------------------------------------------------------------------
% Input data for the problem
%--------------------------------------------------------------------------
% Input files for mesh and FE matrices

sizemesh = load('sizemesh.txt');
sizemesh = sizemesh(end);
config = str2func([mesh.file '_config']);
[flag, param] = config();

timing = containers.Map;

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
        start = cputime;
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
        timing(['FE_', param.study, '_', param.path1]) = cputime - start;
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
                    start = cputime;
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
                    timing(['MDWCAWE_', param.study, '_', param.path2]) = cputime - start;
                end

                %-----------------------------------------------------------------------
                %WCAWE Basis Romain
                %-----------------------------------------------------------------------
                if flag.calculateWCAWE
                    start = cputime;
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
                    timing(['WCAWE_', param.study, '_', param.path2]) = cputime - start;
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
    run_converge(param);
end

if flag.plotMQP
    run_MQP();
end


%--------------------------------------------------------------------------
% convert
%--------------------------------------------------------------------------
if flag.convert2VTK
    run_convert(param, mesh, flag);
end


if flag.calculateTL
    timing = run_STL(param, mesh, timing, sizemesh);
end


if flag.eigen
    eigs(LHS{1},LHS{2},10,'sm')
end



if flag.normalized_error
    normalize_error();
end


timing_PP(timing, param, 'save');


if flag.show_timing
    timing_PP(timing, param, 'display');
end










