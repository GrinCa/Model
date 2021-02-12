function [flag, param] = Modelv5_config()

param.filename = 'Modelv4';

% Input parameters for Matlab calculation
flag.rerun = 1; % to recalculate FreeFem++ matrices
flag.recalculated = 1; % allow WCAWE and/or FE recalculation. if 0 then the
                       % next three flags won't be considered.
flag.calculateFE = 1;  % calculate FE solution, 
flag.calculateMDWCAWE = 0; % calculate MDWCAWE solution
flag.calculateWCAWE = 0; % calculate WCAWE solution


flag.convert2VTK = 0; % convert SOLFE.mat into a .vkt file

flag.eigen = 0;

flag.converge = 0; % this is to post process convergence test, not to perform
                   % one
flag.plotMQP = 0;
flag.calculateTL = 0;
flag.converge_sizemesh = 0;
flag.compare_FE_WCAWE = 0;
flag.normalized_error = 0;

flag.show_timing = 0;


flag.getmatrices = 1; % if 0, the programm won't read matrices, it is really 
% useful if you just post-process the results. Indeed, depending on the
% model, it can be quite long to read the .txt file that contains the
% matrices

if flag.converge || flag.plotMQP || flag.convert2VTK || flag.calculateTL || flag.converge_sizemesh || flag.normalized_error
    flag.getmatrices = 0;
    flag.rerun = 0;
    flag.recalculated = 0;
end





% Material parameters
param.c0 = 340; %m/s
param.rho = 1200; %kg/m3
param.rho0 = 1.29; %kg/m3
%%%%% Background pressure field %%%%%

% Frequency range
param.fmin = 100;
param.fmax = 100;
param.f_range = [param.fmin param.fmax];
param.freqincr = 5; % 20
param.freq = param.fmin : param.freqincr : param.fmax; % frequency range
param.nfreq = length(param.freq);

% Angle range
param.thetamin = 0;
param.thetamax = 0;
param.theta_range = [param.thetamin param.thetamax];
param.thetaincr = 0.05;
param.theta = param.thetamin : param.thetaincr : param.thetamax; % frequency range
param.ntheta = length(param.theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% those frequencies are the frequencies point for Padé expension
param.freqref = [225];
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
param.nvecfreqmin = 3;
param.nvecfreqmax = 3;
param.incrvecfreq = 20;
param.vecfreqrange = param.nvecfreqmin : param.incrvecfreq : param.nvecfreqmax;

param.nvecthetamin = 1;
param.nvecthetamax = 1;
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
P0 = 1;
param.direction = [1;0;1];
param.direction_coeff = param.direction(find(param.direction));

% Definition of the matrices coeffs
param.coeff_LHS = {@(f,theta) 1,@(f,theta) -(2*pi*f)^2};
param.coeff_RHS = @(f,theta,x1,x2) P0*exp((1i*2*pi*f/param.c0).*(param.direction_coeff(1)*x1.*cos(theta)+...
                                                                 param.direction_coeff(2)*x2.*sin(theta)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.matrix_names = ["Kr.txt","Ki.txt","M.txt",...
                      "Hpmlr.txt","Hpmli.txt",...
                      "Qpmlr.txt","Qpmli.txt",...
                      "C1.txt","C2.txt"];
                  
param.study = 'PML';
param.VTK.config = 'quadratic';

end
