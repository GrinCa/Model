function data = post_process(FEmatrices,param,arg,SOLUTION)

    % arg: struct of different argument to perform specific tasks
    % SOLUTION: struct of SOLUTION calculated with FE, MDWCAWE or WCAWE algo
    data = false;

    if strcmp(arg.type,'calculateTL')
        data = calculateTL(FEmatrices,param,SOLUTION);
    end
    
    if strcmp(arg.type,'rel_error')
        data = get_relative_error(arg);
    end
    
    if strcmp(arg.type,'total_pressure')
        data = getTotalPressure(FEmatrices,SOLUTION,param);
    end
    
    if strcmp(arg.type,'normalize_error')
        data = normalize_error(FEmatrices,arg,param);
    end
    
    if strcmp(arg.type,'show_image')
        show_image(arg.matrix);
    end
    
    
end


function p_t = getTotalPressure(FEmatrices,SOL,param)

    ndof = size(FEmatrices.Nodes,1);
    p_t = zeros(ndof,param.nfreq,param.ntheta);
    Scattered_pressure = zeros(ndof,param.nfreq,param.ntheta);
    Scattered_pressure(FEmatrices.BG_PML_nodes,:,:) = SOL(FEmatrices.indexP_BG_PML,:,:);
    for ii=1:param.nfreq
        for jj=1:param.ntheta
            p_t(:,ii,jj) = FEmatrices.BG_pressure(:,ii,jj) + Scattered_pressure(:,ii,jj);
        end
    end
    
end

function TL = calculateTL(FEmatrices,param,SOL)

TL = get_STL(FEmatrices, param, SOL);



end


function Power = calculatePower(surface_matrix,normal_displacement_surface,pressure_surface,param)
    Power = zeros(param.nfreq,param.ntheta);
    for ii=1:param.nfreq
        for jj=1:param.ntheta
            Powertmp = abs(surface_matrix)*pressure_surface(:,ii,jj);
            Vn = 2*pi*param.freq(ii)*normal_displacement_surface(:,ii,jj);
            Power(ii,jj) = 0.5*abs(real(Vn'*Powertmp));
        end
    end
end

function rel_error = get_relative_error(arg)
    rel_error = cell(1,length(arg.APPROX_SOLUTION));
    for ii=1:length(arg.APPROX_SOLUTION)
        rel_error{ii} = log10(abs((arg.APPROX_SOLUTION{ii}-arg.REF_SOLUTION{1})./arg.REF_SOLUTION{1}));
    end
end



function normalized_error = normalize_error(FEmatrices,arg,param)
    % APPROX_SOLUTION is the FE solution whose size = (ndof,nfreq,ntheta)
    % REF_SOLUTION is the FPS solution, same size.
    % normalized_error gives the norm of the relative_error vector calculated
    % with the approx and ref solutions.
    ndof = size(arg.REF_SOLUTION,1);
    normalized_error = zeros(param.nfreq,param.ntheta);
    for ii=1:param.nfreq
       for jj=1:param.ntheta
           rel_error_tmp = (real(arg.APPROX_SOLUTION{1}(:,ii,jj))-real(arg.REF_SOLUTION(:,ii,jj)));%./real(arg.REF_SOLUTION(:,ii,jj));
           normalized_error(ii,jj) = norm(rel_error_tmp)/sqrt(ndof);
       end
    end
end



function show_image(sparse_matrix)

[rows,colums,values] = find(sparse_matrix);
values = ones(length(values),1);

image(real(sparse(rows,colums,values)),'CDataMapping','scaled');

end

   





