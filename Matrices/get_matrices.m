function [FEmatrices,ndof,flag] = get_matrices(flag,mesh,param)

path_soft = getpath();

for ii=1:length(param.matrix_names)
    param.matrix_names{ii} = strcat('Matrices/',mesh.file,'/',param.matrix_names{ii});
end


%--------------------------------------------------------------------------
%Run FreeFem++ script, IF NEEDED
%--------------------------------------------------------------------------

if (flag.rerun) % EDP updated and/or mesh updated
   t_0 = cputime;
   disp('************************');
   disp('*Rerun FreeFem++ script*');
   disp('************************');
   edpcommand = strcat([path_soft.freefem,' ','EDP/', mesh.file, '/', mesh.file,'.edp']);
   system(edpcommand);
   timing.freefem = cputime-t_0;
   disp('*********************************************************');
   output = sprintf('[Get_matrices:infos] CPUtime for building of matrices %.4f s',timing.freefem);
   disp(output);
   disp('*********************************************************');
   flag.recalculated = 1;
end % if

%--------------------------------------------------------------------------
% Get matrices from the files
%--------------------------------------------------------------------------

listLHS = cell(1,length(param.matrix_names)); 


% Matrices of the FE problem
for ii=1:length(param.matrix_names)
    nrows = 0;
    ncolums = 0;
    fid = fopen(param.matrix_names(ii),'rt');
    for jj=1:3
        line = fgets(fid);
        if jj==3
           values = str2num(strtrim(line));
           nrows = values(1);
           ncolums = values(2);
        end
    end
    fclose(fid);
    matrix_data = importdata(param.matrix_names(ii)," ",3);
    matrix_data = matrix_data.data;
    listLHS{ii} = sparse([matrix_data(:,1);nrows-1]+1,[matrix_data(:,2);ncolums-1]+1,[matrix_data(:,3);0]);
end

% Nodes
Nodes = load(strcat(mesh.file,'/',"Nodes.txt"));
ndof = size(Nodes,1);

%--------------------------------------------------------------------------
% return
%--------------------------------------------------------------------------

FEmatrices.Nodes = Nodes; 
FEmatrices = build_global(FEmatrices,listLHS,param,mesh.file);

% RHS
% RHSdata = importdata(strcat(mesh.file,'/',"RHS.txt")," ",3);
% RHSdata = RHSdata.data;
% RHS = diag(sparse(RHSdata(:,1)+1,RHSdata(:,2)+1,RHSdata(:,3)));
% FEmatrices.RHS = RHS(FEmatrices.field);


end


function FEmatrices = build_global(FEmatrices,listLHS,param,FILENAME)
% this function call the specific function for each given problem. You have
% to create a .m file for eash problem called Problem_pattern. For instance
% if Problem1 is one problem to study

pattern = str2func([FILENAME '_pattern']);
FEmatrices = pattern(FEmatrices,listLHS,param,FILENAME);
end


% this following function build the matrices of the system for the plate
% case only. For each new case, we have to build a new function, a pattern
% for partitionning the system.



function FEmatrices = build_BGField(FEmatrices,param)
ndof = size(FEmatrices.Nodes,1);

% Cell of RHS
FEmatrices.RHS_BG = zeros(FEmatrices.size_system,param.nfreq,param.ntheta);

%BG_nodes = 1:1:FEmatrices.size_system;
Field_nodes = FEmatrices.field;

xbg = FEmatrices.Nodes(Field_nodes,1);
ybg = FEmatrices.Nodes(Field_nodes,2);

FEmatrices.BG_pressure = zeros(ndof,param.nfreq,param.ntheta);

P0 = 1;

for ii=1:param.nfreq
    for jj=1:param.ntheta
        k = 2*pi*param.freq(ii)/param.c0;
        BG_Pressure_tmp = P0*exp(1i*k*(xbg*cos(param.theta(jj))+ybg*sin(param.theta(jj))));%propagation (+x,+y), convention exp(-1i*k*x)
        %FEmatrices.BG_pressure(Field_nodes,ii,jj) = BG_Pressure_tmp;
        %Z = FEmatrices.Hbg - (2*pi*param.freq(ii)/param.c0)^2*FEmatrices.Qbg;
        U_inc = zeros(FEmatrices.size_system,1);
        U_inc(FEmatrices.indexfield) = BG_Pressure_tmp;%U_inc(FEmatrices.indexfield) = -Z*BG_Pressure_tmp;
        %U_inc(FEmatrices.indexP_BGPML) = 0;
        FEmatrices.RHS_BG(:,ii,jj) = U_inc;
    end
end
end

































