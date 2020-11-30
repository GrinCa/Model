function [P_inc, P_rad] = get_powers(FEmatrices, param, SOL)

[rayleigh_matrices, element_data, normal_direction,Plan2D] = RIM(FEmatrices, param);
P_rad = zeros(param.nfreq, param.ntheta);
P_inc = zeros(param.nfreq, param.ntheta);

% Uortho: normal displacement of the plate
Uortho = zeros(size(FEmatrices.Nodes,1),param.nfreq,param.ntheta);
Uortho(FEmatrices.plate_nodes,:,:) = SOL(normal_direction:3:3*length(FEmatrices.plate_nodes),:,:);

Uortho_element = zeros(size(element_data,1),param.nfreq,param.ntheta);
for kk=1:length(element_data)
    for ii=1:param.nfreq
        for jj=1:param.ntheta
            elements = element_data(kk,:);
            Uortho_element(kk,ii,jj) = mean(Uortho(elements,ii,jj));
        end
    end
end

% Calculation of P_rad
for ii=1:param.nfreq
    for jj=1:param.ntheta
        P_rad(ii,jj) = (2*pi*param.freq(ii))^2*Uortho_element(:,ii,jj)'*rayleigh_matrices{ii}*Uortho_element(:,ii,jj);
    end
end

for ii=1:param.nfreq
    for jj=1:param.ntheta
        Pressure_surface = FEmatrices.RHScoeffderiv_fun{1,1,1}(param.freq(ii),...
                                                               param.theta(jj),...
                                                               FEmatrices.Nodes(FEmatrices.PlateIn,Plan2D(1)),...
                                                               FEmatrices.Nodes(FEmatrices.PlateIn,Plan2D(2)));
        P_inc(ii,jj) = (2*pi*param.freq(ii))^2*Uortho(FEmatrices.PlateIn,ii,jj)'*...
                        FEmatrices.Surf_matrix*Pressure_surface;
    end
end

end