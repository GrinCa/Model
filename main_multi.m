%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Test case using WCAWE                             %
%                                                                         %
%                             March 2020                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Init main program
%--------------------------------------------------------------------------

clear FEmatrices param flag;
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
flag.rerun = 0; % to recalculate FreeFem++ matrices
flag.recalculated = 1; % allow WCAWE and/or FE recalculation
flag.calculateFE = 0;  % calculate FE solution
flag.calculateMDWCAWE = 1; % calculate MDWCAWE solution
flag.calculateWCAWE = 0; % calculate WCAWE solution

flag.plotcomparison = 0; % plot comparison between FE and WCAWE
flag.comparisonMULTI = 0;
flag.show_graph = 0;

flag.converge = 0;
flag.convert2VTK = 1; % convert SOLFE.mat into a .vkt file
flag.plotMQP = 0;

flag.getmatrices = 1;

if flag.converge || flag.plotMQP || flag.convert2VTK
    flag.getmatrices = 0;
    flag.rerun = 0;
    flag.recalculated = 0;
end



% Input files for mesh and FE matrices
mesh.file = 'Debug';
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
param.fmin = 50;
param.fmax = 150;
param.f_range = [param.fmin param.fmax];
param.freqincr = 10; % 20
param.freq = param.fmin : param.freqincr : param.fmax; % frequency range
param.nfreq = length(param.freq);

% Angle range
param.thetamin = -1;
param.thetamax = 1;
param.theta_range = [param.thetamin param.thetamax];
param.thetaincr = 0.2;
param.theta = param.thetamin : param.thetaincr : param.thetamax; % frequency range
param.ntheta = length(param.theta);

P0 = 1;
param.direction = [1;0;1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% those frequencies are the frequencies point for Padé expension
param.freqref = [70 105 130];
param.nfreqref = length(param.freqref);

param.thetaref = [0];
param.nthetaref = length(param.thetaref);

% Input data for the loop over expansion orders. Note that for each
% frequency sweep the number of vectors from the WCAWE basis will depend on
% the number of point for Padé expension. For instance, if we have 2 points
% for expansion, and nvecfreq=5 (order of expansion), we will have 15
% vectors in the basis, 5 by intervals.
param.nvecfreqmin = 3;
param.nvecfreqmax = 3;
param.incrvec = 20;
param.vecfreqrange = param.nvecfreqmin : param.incrvec : param.nvecfreqmax;


param.nvecthetamin = 4;
param.nvecthetamax = 4;
param.vecthetarange = param.nvecthetamin : param.incrvec : param.nvecthetamax;

%Identificator
folder.idData = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']['...
                num2str(int16(180*param.theta_range(1)/pi)) '_' num2str(int16(180*param.theta_range(2)/pi)) ']/'...
                '[' num2str(param.nvecfreqmin) '_' num2str(param.nvecthetamin) ']['...
                replace(num2str(param.freqref),'  ','_') '_' replace(num2str(param.thetaref),'  ','_') ']'];
folder.splitData = 16;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the matrices coeffiscients
coeff_LHS = {@(f,theta) 1,@(f,theta) -(2*pi*f)^2};
coeff_RHS = @(f,theta,x1,x2) P0/800*exp(1i*(2*pi*f/param.c0).*(x1.*cos(theta)-x2.*sin(theta)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generation of the different folder to store data if they don't already exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genfolders(mesh,folder);
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
                    "H.txt","Q.txt",...
                    "C.txt"];

    [FEmatrices,ndof,timing,flag] = get_matrices(timing,flag,mesh,matrix_names,param);
    Nodes = FEmatrices.Nodes;
    LHS = FEmatrices.LHS;
    nLHS = length(LHS);
    param.idx_out = 1:1:FEmatrices.size_system;
    % Save data only for FE solution
    save(['Matrices/',mesh.file,'/',folder.idData(1:folder.splitData),'/','DATA_sizemesh_',num2str(sizemesh),'.mat'],'FEmatrices','param','timing');
end





if flag.recalculated
    
% Initialize array of RHSderiv, Cell array linear comb coefficients derivative functions
deriv_deg = [param.nvecfreqmax, param.nvecthetamax];

if exist('Derivatives/derivative_orders.mat','file') ~= 2
    disp('*************************************');
    disp('* Recalculate all cross derivatives *');
    disp('*************************************');
    create_cross_derivatives(LHS,coeff_LHS,...
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

% get function handle of cross derivatives
[coeff_deriv_fun,RHScoeffderiv_fun] = get_coeff_deriv_matrices([param.nvecfreqmax,param.nvecthetamax],nLHS);
 
%--------------------------------------------------------------------------
% Calculate reference finite element parametric sweep
%--------------------------------------------------------------------------


if flag.calculateFE == 1
   disp('**************************');
   disp('* Compute FE calculation *');
   disp('**************************');
   t_FE = cputime;
   SOLFE = zeros(FEmatrices.size_system,param.nfreq,param.ntheta); %size ndof,nfreq
   % Parametric loop calculation
   for ii=1:param.nfreq
       for jj=1:param.ntheta
           tic;
           disp(['[FE] [Frequency, Theta] = [',num2str(param.freq(ii)),',',num2str(param.theta(jj)/pi*180),']']);
           Aglob = sparse(size(LHS{1},1),size(LHS{1},2));
           for kk = 1:nLHS
              Aglob = Aglob + coeff_LHS{kk}(param.freq(ii))*LHS{kk};
           end
           RHS_0 = build_RHS(param.freq(ii),param.theta(jj),FEmatrices,RHScoeffderiv_fun,[1,1],param.direction);
           FEmatrices.RHS_0 = RHS_0;
           resp_P = Aglob\RHS_0; %FEmatrices.RHS;
           SOLFE(:,ii,jj) = resp_P;
           toc;
       end
   end
   timing.FE = cputime-t_FE;
end



%--------------------------------------------------------------------------
% Calculate WCAWE parametric sweep
%--------------------------------------------------------------------------

if flag.calculateMDWCAWE || flag.calculateWCAWE
    
for nvecfreq=param.vecfreqrange
    for nvectheta=param.vecthetarange
        %%%%%
        param.nvecfreq = nvecfreq;
        param.nvectheta = nvectheta;
        %%%%%
        [LHScoeffderiv_fun,RHScoeffderiv_fun] = get_coeff_deriv_matrices([nvecfreq,nvectheta],nLHS);
        FEmatrices.LHScoeffderiv_fun = LHScoeffderiv_fun;
        FEmatrices.RHScoeffderiv_fun = RHScoeffderiv_fun;
        [LHScoeffderiv,RHSderiv] = fill_array_WCAWE(FEmatrices,LHScoeffderiv_fun,RHScoeffderiv_fun,param);
        
        if flag.calculateMDWCAWE
        %-----------------------------------------------------------------------
        %MDWCAWE Basis Julien
        %-----------------------------------------------------------------------
        disp('****************************************');
        disp('* Compute calculation of MDWCAWE basis *');
        disp('****************************************');
        t_MDWCAWE = cputime;
        MDWCAWE = [];
        for n=1:param.nfreqref
            for m=1:param.nthetaref
                MDWCAWEtmp = MDWCAWE_basis(LHS,LHScoeffderiv{n,m},RHSderiv{n,m},nvecfreq,nvectheta);
                MDWCAWE = [MDWCAWE MDWCAWEtmp];
            end
        end
        [uu,vv,ww] = svd(MDWCAWE,0);
        iiselect = find(diag(vv)>vv(1,1)*1e-15);
        MDWCAWEsvd = uu(:,iiselect);
        nsvd = size(MDWCAWEsvd,2);
        output = sprintf("[SVD:INFO [MDWCAWE]] Number of selected vectors %d/%d\n",nsvd,size(MDWCAWE,2));
        disp(output);
        disp('**** Compute MDWCAWE projection ****');
        SOLMDWCAWE = Solve_MDWCAWE(FEmatrices,MDWCAWEsvd,param);
        disp('**** MDWCAWE projection done ****');
        timing.MDWCAWE = cputime-t_MDWCAWE;
        end
        %-----------------------------------------------------------------------
        %WCAWE Basis Romain
        %-----------------------------------------------------------------------
        if flag.calculateWCAWE
            disp('**************************************');
            disp('* Compute calculation of WCAWE basis *');
            disp('**************************************');
            t_WCAWE = cputime;
            % Setup list of derivative coefficients and RHS derivatives
            WCAWE = [];
            Wtranstmp_order = param.nvecfreq+param.nvectheta-1;
            for n=1:param.nfreqref
                for m=1:param.nthetaref
                    LHScoeffderivtmp = zeros(nLHS,Wtranstmp_order);
                    RHSderivtmp = zeros(FEmatrices.size_system,Wtranstmp_order);
                    counter = 1;
                    for ii = 1:nvecfreq
                        if ii<nvecfreq
                            LHScoeffderivtmp(:,counter) = LHScoeffderiv{n,m}(:,ii,1);
                            RHSderivtmp(:,counter) = RHSderiv{n,m}(:,ii,1);
                            counter = counter+1;
                        elseif ii==nvecfreq
                            for jj = 1:nvectheta
                               LHScoeffderivtmp(:,counter) = LHScoeffderiv{n,m}(:,ii,jj);
                               RHSderivtmp(:,counter) = RHSderiv{n,m}(:,ii,jj);
                               counter = counter+1;
                            end % jj
                        end
                    end % ii
                    WCAWEtmp = WCAWE_basis(FEmatrices,LHScoeffderivtmp,RHSderivtmp,Wtranstmp_order);
                    WCAWE = [WCAWE WCAWEtmp];
                    
                    LHScoeffderivtmp = zeros(nLHS,Wtranstmp_order);
                    RHSderivtmp = zeros(FEmatrices.size_system,Wtranstmp_order);
                    counter = 1;
                    for jj = 1:nvectheta
                        if jj<nvectheta
                            LHScoeffderivtmp(:,counter) = LHScoeffderiv{n,m}(:,1,jj);
                            RHSderivtmp(:,counter) = RHSderiv{n,m}(:,1,jj);
                            counter = counter+1;
                        elseif jj==nvectheta
                            for ii = 1:nvecfreq
                               LHScoeffderivtmp(:,counter) = LHScoeffderiv{n,m}(:,ii,jj);
                               RHSderivtmp(:,counter) = RHSderiv{n,m}(:,ii,jj);
                               counter = counter+1;
                            end % ii
                        end
                    end % jj
                    WCAWEtmp = WCAWE_basis(FEmatrices,LHScoeffderivtmp,RHSderivtmp,Wtranstmp_order);
                    WCAWE = [WCAWE WCAWEtmp];
                end
            end
            [uu,vv,ww] = svd(WCAWE,0);
            iiselect = find(diag(vv)>vv(1,1)*1e-15);
            WCAWEsvd = uu(:,iiselect);
            nsvd = size(WCAWEsvd,2);
            output = sprintf("[SVD:INFO [WCAWE]] Number of selected vectors %d/%d\n",nsvd,size(MDWCAWE,2));
            disp(output);
            disp('**** Compute WCAWE projection ****');
            SOLWCAWE = Solve_MDWCAWE(FEmatrices,WCAWEsvd,param);
            disp('**** WCAWE projection done ****');
            timing.WCAWE = cputime-t_WCAWE;
        end
        
        %--------------------------------------------------------------------------
        % Saves
        %--------------------------------------------------------------------------
        % FE solution
        if flag.calculateFE
            save(['Matrices/',mesh.file,'/',folder.idData(1:folder.splitData),'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat'],'SOLFE');
        end
        if flag.calculateMDWCAWE
            id_sample = ['_sizemesh_',num2str(sizemesh)];
            save(['Matrices/',mesh.file,'/',folder.idData,'/SOLMDWCAWE',id_sample,'.mat'],'SOLMDWCAWE');
        end
        if flag.calculateWCAWE
            id_sample = ['_sizemesh_',num2str(sizemesh)];
            save(['Matrices/',mesh.file,'/',folder.idData,'/SOLWCAWE',id_sample,'.mat'],'SOLWCAWE');
        end
        % save FEmatrices which contains all the data of the simulation for
        % each mesh
        save(['Matrices/',mesh.file,'/',folder.idData,'/','DATA_sizemesh_',num2str(sizemesh),'.mat'],'FEmatrices','param','timing');

    end
end
end
end

%--------------------------------------------------------------------------
% Post processing
%--------------------------------------------------------------------------

if flag.show_graph
    show_graph('compare_results',mesh,param,folder);
end

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
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',folder.idData(1:folder.splitData),'/','DATA_sizemesh_',num2str(sizemesh),'.mat']));
    FEmatrices = DATA{1};
    param = DATA{2};
    %%%
    PARTITION = cell(2);
    
    PARTITION{1} = {'scalar',...
                    FEmatrices.indexp,...
                    FEmatrices.cavity_nodes,...
                    'CAVITY_PRESSURE'};
    PARTITION{2} = {'vector',...
                    [FEmatrices.indexu1;...
                     FEmatrices.indexu2;...
                     FEmatrices.indexu3],...
                     FEmatrices.plate_nodes,...
                     'DISPLACEMENT'};
    %%%
    range = {6:1:6, 6:1:6};
    
    if flag.calculateMDWCAWE
        SOLMDWCAWE = struct2cell(load(['Matrices/',mesh.file,'/',folder.idData,'/SOLMDWCAWE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLMDWCAWE = SOLMDWCAWE{1};
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLMDWCAWE,PARTITION,param,range);
    end
    if flag.calculateWCAWE
        SOLWCAWE = struct2cell(load(['Matrices/',mesh.file,'/',folder.idData,'/SOLWCAWE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLWCAWE = SOLWCAWE{1};
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLWCAWE,PARTITION,param,range);
    end
    if flag.calculateFE
        SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',folder.idData(1:folder.splitData),'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLFE = SOLFE{1};
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLFE,PARTITION,param,range);
    end
end














