function [P_inc, P_rad] = get_powers(FEmatrices, param, SOL)

%--------------------------------------------------------------------------
% Calulation of P_inc

% Uortho: normal displacement of the plate
Uortho = zeros(size(FEmatrices.Nodes,1),param.nfreq,param.ntheta);
Uortho(FEmatrices.plate_nodes,:,:) = SOL(FEmatrices.normal_direction:3:3*length(FEmatrices.plate_nodes),:,:);
P_inc = zeros(param.nfreq, param.ntheta);
for ii=1:param.nfreq
    for jj=1:param.ntheta
        Pressure_surface = FEmatrices.RHScoeffderiv_fun{1,1,1}(param.freq(ii),...
                                                               param.theta(jj),...
                                                               FEmatrices.Nodes(FEmatrices.PlateIn,FEmatrices.plan2D(1)),...
                                                               FEmatrices.Nodes(FEmatrices.PlateIn,FEmatrices.plan2D(2)));
        P_inc(ii,jj) = 1i*2*pi*param.freq(ii)*Uortho(FEmatrices.PlateIn,ii,jj)'*...
                       FEmatrices.SurfIn_matrix*Pressure_surface;
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Calculation of P_rad
if strcmp(param.STL.config, 'PML')
    P_rad = PML_method(FEmatrices, param, SOL);
elseif strcmp(param.STL.config, 'RAYLEIGH')
    P_rad = RAYLEIGH_method(FEmatrices, param, SOL);
end
%--------------------------------------------------------------------------


end

function P_rad = PML_method(FEmatrices, param, SOL)

P_rad = zeros(param.nfreq, param.ntheta);
% Uortho: normal displacement of the plate
Uortho = zeros(size(FEmatrices.Nodes,1),param.nfreq,param.ntheta);
Uortho(FEmatrices.plate_nodes,:,:) = SOL(FEmatrices.normal_direction:3:3*length(FEmatrices.plate_nodes),:,:);

pressure = SOL(FEmatrices.indexPlateExt,:,:);

for ii=1:param.nfreq
    for jj=1:param.ntheta
        P_rad(ii,jj) = 1i*2*pi*param.freq(ii)*Uortho(FEmatrices.PlateExt,ii,jj)'*...
                            FEmatrices.SurfExt_matrix*pressure(:,ii,jj); % v = -iwUortho and 
                        % v^H = conjugate(v^T) , therefore:
                        % v^H = iwUortho^H
    end
end

end


function P_rad = RAYLEIGH_method(FEmatrices, param, SOL)

[rayleigh_matrices, element_data] = RIM(FEmatrices, param);
P_rad = zeros(param.nfreq, param.ntheta);


% Uortho: normal displacement of the plate
Uortho = zeros(size(FEmatrices.Nodes,1),param.nfreq,param.ntheta);
Uortho(FEmatrices.plate_nodes,:,:) = SOL(FEmatrices.normal_direction:3:3*length(FEmatrices.plate_nodes),:,:);

% plot3(FEmatrices.Nodes(FEmatrices.plate_nodes,FEmatrices.plan2D(1)),...
%       FEmatrices.Nodes(FEmatrices.plate_nodes,FEmatrices.plan2D(2)),...
%       real(Uortho(FEmatrices.plate_nodes,30,1)),'+');

Uortho_element = zeros(size(element_data,1),param.nfreq,param.ntheta);
for kk=1:length(element_data)
    for ii=1:param.nfreq
        for jj=1:param.ntheta
            Uortho_element(kk,ii,jj) = mean(Uortho(element_data(kk,:),ii,jj));
        end
    end
end

% Calculation of P_rad
for ii=1:param.nfreq
    for jj=1:param.ntheta
        P_rad(ii,jj) = -(2*pi*param.freq(ii))^2*Uortho_element(:,ii,jj)'*rayleigh_matrices{ii}*Uortho_element(:,ii,jj);
    end
end



end