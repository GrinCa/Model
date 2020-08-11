%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Test case using WCAWE                             %
%                                                                         %
%                             March 2020                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Init main program
%--------------------------------------------------------------------------

%clear FEmatrices param flag;
warning('off', 'MATLAB:nearlySingularMatrix');

%--------------------------------------------------------------------------
% Folders
%--------------------------------------------------------------------------

% Add folders for Mumps, WCAWE, and Mesh functions/files
global Meshfold WCAWEfold Mumps DataMap Derivatives
Meshfold = 'Matrices';
WCAWEfold = 'WCAWE';
Mumps = 'Mumps';
DataMap = 'DataMap';
Derivatives = 'Derivatives';



addpath(genpath(strcat(pwd,'/',Meshfold)));
addpath(genpath(strcat(pwd,'/',WCAWEfold)));
addpath(genpath(strcat(pwd,'/',Mumps)));
addpath(genpath(strcat(pwd,'/',DataMap)));
addpath(genpath(strcat(pwd,'/',Derivatives)));




%--------------------------------------------------------------------------
% Input data for the problem
%--------------------------------------------------------------------------

% Input parameters for Matlab calculation
flag.rerun = 1; % to recalculate FreeFem++ matrices
flag.recalculated = 1; % allow WCAWE and/or FE recalculation
flag.calculateFE = 1;  % calculate FE solution
flag.calculateMDWCAWE = 0; % calculate MDWCAWE solution
flag.calculateWCAWE = 0; % calculate WCAWE solution



flag.converge = 0;
flag.convert2VTK = 0; % convert SOLFE.mat into a .vkt file
flag.plotMQP = 0;
flag.calculateTL = 0;
flag.converge_sizemesh = 0;
flag.compare_FE_WCAWE = 0;
flag.normalized_error = 0;

flag.getmatrices = 1;

if flag.converge || flag.plotMQP || flag.convert2VTK || flag.calculateTL || flag.converge_sizemesh || flag.normalized_error
    flag.getmatrices = 0;
    flag.rerun = 0;
    flag.recalculated = 0;
end



% Input files for mesh and FE matrices
mesh.file = 'Platev2';
sizemesh = load('sizemesh.txt');
sizemesh = sizemesh(end);

% define timing struct to access to time calculation of each method                                                    
timing.freefem = 0;
timing.WCAWE = 0;
timing.MDWCAWE = 0;
timing.FE = 0;                                                    


% Material parameters
param.rho = 1.2;
param.rhoS = 1200;
param.c0 = 340;

%%%%% Background pressure field %%%%%

% Frequency range
param.fmin = 600;
param.fmax = 600;
param.f_range = [param.fmin param.fmax];
param.freqincr = 2; % 20
param.freq = param.fmin : param.freqincr : param.fmax; % frequency range
param.nfreq = length(param.freq);

% Angle range
param.thetamin = 0;
param.thetamax = 0;
param.theta_range = [param.thetamin param.thetamax];
param.thetaincr = 0.05;
param.theta = param.thetamin : param.thetaincr : param.thetamax; % frequency range
param.ntheta = length(param.theta);

P0 = 1;
param.direction = [1;1;0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% those frequencies are the frequencies point for Padé expension
param.freqref = [400];
param.nfreqref = length(param.freqref);

param.thetaref = [0];
param.nthetaref = length(param.thetaref);

% interval_construct enables us to build sub basis for WCAWE by
% by using the ref frequencies as we want. We can choose which ref freq to
% add for each sub basis. The number of sub basis for WCAWE before SVD is
% equal to length(interval_construct)
param.interval_construct = {{[1]};
                            {[1]}};

% Input data for the loop over expansion orders. Note that for each
% frequency sweep the number of vectors from the WCAWE basis will depend on
% the number of point for Padé expension. For instance, if we have 2 points
% for expansion, and nvecfreq=5 (order of expansion), we will have 15
% vectors in the basis, 5 by intervals.
param.nvecfreqmin = 30;
param.nvecfreqmax = 30;
param.incrvecfreq = 20;
param.vecfreqrange = param.nvecfreqmin : param.incrvecfreq : param.nvecfreqmax;

param.nvecthetamin = 10;
param.nvecthetamax = 10;
param.incrvectheta = 20;
param.vecthetarange = param.nvecthetamin : param.incrvectheta : param.nvecthetamax;

%Identificator
param.path1 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']['...
                num2str(int16(180*param.theta_range(1)/pi)) '_' num2str(int16(180*param.theta_range(2)/pi)) ']'];
param.path2 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']['...
                num2str(int16(180*param.theta_range(1)/pi)) '_' num2str(int16(180*param.theta_range(2)/pi)) ']/'...
                '[' num2str(param.nvecfreqmin) '_' num2str(param.nvecthetamin) ']['...
                replace(num2str(param.freqref),' ','_') '][' replace(num2str(param.thetaref),' ','_') ']'];
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the matrices coeffiscients
coeff_LHS = {@(f,theta) 1,@(f,theta) -(2*pi*f)^2};
coeff_RHS = @(f,theta,x1,x2) P0*exp(1i*(2*pi*f/param.c0).*(x1.*cos(theta)+x2.*sin(theta)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generation of the different folder to store data if they don't already exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param = genfolders(mesh,param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% Build intervals for the Parametric sweep
%--------------------------------------------------------------------------
param = build_interval(param);

%--------------------------------------------------------------------------
% Matrices calculated with Freefem++
%--------------------------------------------------------------------------
if flag.getmatrices
    matrix_names = ["K.txt","M.txt",... % matrices defined on the incident domain
                    "Hr.txt","Hi.txt",...
                    "Qr.txt","Qi.txt",...
                    "C1.txt","C2.txt"];
                
%     matrix_names = ["K.txt","M.txt",... % matrices defined on the incident domain
%                     "Hbg.txt","Qbg.txt",...
%                     "Hcavr.txt","Qcavr.txt",...
%                     "Hcavi.txt","Qcavi.txt",...
%                     "Hpmlr.txt","Hpmli.txt",...
%                     "Qpmlr.txt","Qpmli.txt",...
%                     "C1.txt","C2.txt"];

                
    [FEmatrices,ndof,timing,flag] = get_matrices(timing,flag,mesh,matrix_names,param);
    Nodes = FEmatrices.Nodes;
    LHS = FEmatrices.LHS;
    nLHS = length(LHS);
    param.idx_out = 1:1:FEmatrices.size_system;
    % Save data only for FE solution
    save(['Matrices/',mesh.file,'/',param.path1,'/','DATA_sizemesh_',num2str(sizemesh),'.mat'],'FEmatrices','param','timing');
end



if flag.recalculated

    %--------------------------------------------------------------------------
    % Calculate reference finite element parametric sweep
    %--------------------------------------------------------------------------
    
    [FEmatrices.LHScoeffderiv_fun, FEmatrices.RHScoeffderiv_fun] = get_coeff_deriv_matrices();

    if flag.calculateFE == 1
       disp('**************************');
       disp('* Compute FE calculation *');
       disp('**************************');
       t_FE = cputime;
       SOLFE = zeros(FEmatrices.size_system,param.nfreq,param.ntheta); %size ndof,nfreq
       % Parametric loop calculation
       id = initmumps;
       id = zmumps(id);
       id.JOB = 1; %analysis
       matrix_analysis = LHS{1} + LHS{2};
       id = zmumps(id,matrix_analysis); % perform analysis
       for ii=1:param.nfreq
           for jj=1:param.ntheta
               tic;
               disp(['[FE] [Frequency, Theta] = [',num2str(param.freq(ii)),',',num2str(param.theta(jj)/pi*180),']']);
               Aglob = sparse(size(LHS{1},1),size(LHS{1},2));
               for kk = 1:nLHS
                  Aglob = Aglob + coeff_LHS{kk}(param.freq(ii))*LHS{kk};
               end

               id.JOB = 5;
               id.RHS = build_RHS(param.freq(ii),param.theta(jj),FEmatrices,[1,1],param);
               id = zmumps(id,Aglob);
               resp_P = id.SOL;
               SOLFE(:,ii,jj) = resp_P;
               toc;
           end
       end
       id.JOB = -2; id = zmumps(id); %delete the instance of Mumps
       timing.FE = cputime-t_FE;
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
            create_cross_derivatives(nLHS,coeff_LHS,...
                                     coeff_RHS,deriv_deg,...
                                     {'f';'theta'},{'f,theta';...
                                                    'f,theta';...
                                                    'f,theta,x1,x2'});
        end
        load('Derivatives/derivative_orders.mat');
        if ~isempty(find(derivative_orders-deriv_deg<0))
            disp('*************************************');
            disp('* Recalculate all cross derivatives *');
            disp('*************************************');
            create_cross_derivatives(LHS,coeff_LHS,...
                                     coeff_RHS,deriv_deg,...
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
                [FEmatrices.LHScoeffderiv_fun, FEmatrices.RHScoeffderiv_fun] = get_coeff_deriv_matrices([nvecfreq,nvectheta],nLHS);
                [FEmatrices.LHScoeffderiv, FEmatrices.RHSderiv] = fill_array_WCAWE(FEmatrices,param);

                %-----------------------------------------------------------------------
                %MDWCAWE Basis Julien
                %-----------------------------------------------------------------------
                if flag.calculateMDWCAWE
                    disp('****************************************');
                    disp('* Compute calculation of MDWCAWE basis *');
                    disp('****************************************');
                    t_MDWCAWE = tic;
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
                    timing.MDWCAWE = toc;
                end

                %-----------------------------------------------------------------------
                %WCAWE Basis Romain
                %-----------------------------------------------------------------------
                if flag.calculateWCAWE
                    disp('**************************************');
                    disp('* Compute calculation of WCAWE basis *');
                    disp('**************************************');
                    t_WCAWE = tic;
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
                    timing.WCAWE = toc;
                end
            end
       end
    end
% if flag.calculateMDWCAWE || flag.calculateWCAWE
%     
% for nvecfreq=param.vecfreqrange
%     for nvectheta=param.vecthetarange
%         %%%%%
%         param.nvecfreq = nvecfreq;
%         param.nvectheta = nvectheta;
%         %%%%%
%         [LHScoeffderiv_fun,RHScoeffderiv_fun] = get_coeff_deriv_matrices([nvecfreq,nvectheta],nLHS);
%         FEmatrices.LHScoeffderiv_fun = LHScoeffderiv_fun;
%         FEmatrices.RHScoeffderiv_fun = RHScoeffderiv_fun;
%         [LHScoeffderiv,RHSderiv] = fill_array_WCAWE(FEmatrices,param);
%         
%         if flag.calculateMDWCAWE
%         %-----------------------------------------------------------------------
%         %MDWCAWE Basis Julien
%         %-----------------------------------------------------------------------
%         disp('****************************************');
%         disp('* Compute calculation of MDWCAWE basis *');
%         disp('****************************************');
%         t_MDWCAWE = cputime;
%         MDWCAWE = [];
%         for n=1:param.nfreqref
%             for m=1:param.nthetaref
%                 MDWCAWEtmp = MDWCAWE_basis(LHS,LHScoeffderiv{n,m},RHSderiv{n,m},nvecfreq,nvectheta);
%                 MDWCAWE = [MDWCAWE MDWCAWEtmp];
%             end
%         end
%         [uu,vv,ww] = svd(MDWCAWE,0);
%         iiselect = find(diag(vv)>vv(1,1)*1e-15);
%         MDWCAWEsvd = uu(:,iiselect);
%         nsvd = size(MDWCAWEsvd,2);
%         output = sprintf("[MDWCAWE:INFO][SVD] Number of selected vectors %d/%d\n",nsvd,size(MDWCAWE,2));
%         disp(output);
%         disp('**** Compute MDWCAWE projection ****');
%         SOLMDWCAWE = Solve(FEmatrices,MDWCAWEsvd,param);
%         disp('**** MDWCAWE projection done ****');
%         timing.MDWCAWE = cputime-t_MDWCAWE;
%         end
%         %-----------------------------------------------------------------------
%         %WCAWE Basis Romain
%         %-----------------------------------------------------------------------
%         if flag.calculateWCAWE
%             disp('**************************************');
%             disp('* Compute calculation of WCAWE basis *');
%             disp('**************************************');
%             t_WCAWE = cputime;
%             % Setup list of derivative coefficients and RHS derivatives
%             WCAWE = [];
%             Wtranstmp_order = param.nvecfreq+param.nvectheta-1;
%             for n=1:param.nfreqref
%                 for m=1:param.nthetaref
%                     LHScoeffderivtmp = zeros(nLHS,Wtranstmp_order);
%                     RHSderivtmp = zeros(FEmatrices.size_system,Wtranstmp_order);
%                     counter = 1;
%                     for ii = 1:nvecfreq
%                         if ii<nvecfreq
%                             LHScoeffderivtmp(:,counter) = LHScoeffderiv{n,m}(:,ii,1);
%                             RHSderivtmp(:,counter) = RHSderiv{n,m}(:,ii,1);
%                             counter = counter+1;
%                         elseif ii==nvecfreq
%                             for jj = 1:nvectheta
%                                LHScoeffderivtmp(:,counter) = LHScoeffderiv{n,m}(:,ii,jj);
%                                RHSderivtmp(:,counter) = RHSderiv{n,m}(:,ii,jj);
%                                counter = counter+1;
%                             end % jj
%                         end
%                     end % ii
%                     WCAWEtmp = WCAWE_basis(FEmatrices,LHScoeffderivtmp,RHSderivtmp,Wtranstmp_order);
%                     WCAWE = [WCAWE WCAWEtmp];
%                     
%                     LHScoeffderivtmp = zeros(nLHS,Wtranstmp_order);
%                     RHSderivtmp = zeros(FEmatrices.size_system,Wtranstmp_order);
%                     counter = 1;
%                     for jj = 1:nvectheta
%                         if jj<nvectheta
%                             LHScoeffderivtmp(:,counter) = LHScoeffderiv{n,m}(:,1,jj);
%                             RHSderivtmp(:,counter) = RHSderiv{n,m}(:,1,jj);
%                             counter = counter+1;
%                         elseif jj==nvectheta
%                             for ii = 1:nvecfreq
%                                LHScoeffderivtmp(:,counter) = LHScoeffderiv{n,m}(:,ii,jj);
%                                RHSderivtmp(:,counter) = RHSderiv{n,m}(:,ii,jj);
%                                counter = counter+1;
%                             end % ii
%                         end
%                     end % jj
%                     WCAWEtmp = WCAWE_basis(FEmatrices,LHScoeffderivtmp,RHSderivtmp,Wtranstmp_order);
%                     WCAWE = [WCAWE WCAWEtmp(:,2:end-1)];
%                 end
%             end
%             [uu,vv,ww] = svd(WCAWE,0);
%             iiselect = find(diag(vv)>vv(1,1)*1e-15);
%             WCAWEsvd = uu(:,iiselect);
%             nsvd = size(WCAWEsvd,2);
%             output = sprintf("[WCAWE:INFO][SVD] Number of selected vectors %d/%d\n",nsvd,size(WCAWE,2));
%             disp(output);
%             disp('**** Compute WCAWE projection ****');
%             SOLWCAWE = Solve(FEmatrices,WCAWEsvd,param);
%             disp('**** WCAWE projection done ****');
%             timing.WCAWE = cputime-t_WCAWE;
%         end
%     end
% end
% end

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
    %%%
%     PARTITION = cell(5);
%     
%     PARTITION{1} = {'scalar',...
%                     FEmatrices.indexP_CAVITY,...
%                     FEmatrices.cavity_nodes,...
%                     'CAVITY_PRESSURE'};
%     PARTITION{2} = {'scalar',...
%                      FEmatrices.indexP_BG_PML,...
%                      FEmatrices.BG_PML_nodes,...
%                      'SCATTERED_PRESSURE'};
%     PARTITION{3} = {'vector',...
%                     [FEmatrices.indexu1;...
%                      FEmatrices.indexu2;...
%                      FEmatrices.indexu3],...
%                      FEmatrices.plate_nodes,...
%                      'DISPLACEMENT'};
%     PARTITION{4} = {'data',...
%                      FEmatrices.BG_pressure,...
%                      'BG_PRESSURE'};
%     Force_vector = zeros(size(FEmatrices.Nodes,1),param.nfreq,param.ntheta);
%     Force_vector(FEmatrices.BG_nodes,:,:) = FEmatrices.RHS_BG(FEmatrices.indexP_BG,:,:);
%     PARTITION{5} = {'data',...
%                      Force_vector,...
%                      'FORCE_VECTOR'};

     PARTITION = cell(2);
     PARTITION{1} = {'scalar',...
                     FEmatrices.indexP,...
                     FEmatrices.acoustic_nodes,...
                    'PRESSURE'};
     PARTITION{2} = {'vector',...
                    [FEmatrices.indexu1;...
                     FEmatrices.indexu2;...
                     FEmatrices.indexu3],...
                     FEmatrices.elastic_nodes,...
                     'DISPLACEMENT'};
                              

    %%%
    range = {1:1:1, 1:1:1};

    if flag.calculateFE
        SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLFE = SOLFE{1};
%         arg.type = 'total_pressure';
%         PARTITION{5} = {'data',...
%                         post_process(FEmatrices,param,arg,SOLFE),...
%                         'Total_pressure'};
        
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLFE,PARTITION,param,range);
    end
end

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
   
    argcomp.zlabel = 'TL';
    argcomp.title = 'MDWCAWE';
    argcomp.type = 'plotTL';
    argcomp.split = 0;
    argcomp.name_plot = 'Comparison_FE_MDWCAWE';
    argcomp.label = VALUES_name;
    argcomp.external_plot.is_needed = false;
    show_graph(argcomp,VALUES_SOL(1,:),mesh,param);
    
    argcomp.zlabel = 'TL';
    argcomp.title = 'WCAWE';
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
    argcomp.zlabel = 'relative error';
    argcomp.title = '';
    argcomp.split = 1;
    argcomp.name_plot = 'Relative_error';
    argcomp.label = VALUES_name(2:end);
    argcomp.external_plot.is_needed = false;
    show_graph(argcomp,post_process(FEmatrices,param,argtmp),mesh,param);
    
    
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













