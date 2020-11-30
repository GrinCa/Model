function [LHScoeffderiv,RHSderiv] = fill_array_WCAWE(FEmatrices,param)

% function which automatically fill the array of matrices coeff
coeff_derivgen_fun = @(freq,theta) cellfun(@(cellfunc) cellfunc(freq,theta),FEmatrices.LHScoeffderiv_fun);

% Cell which contain array of coeff (LHS and LHS) for each ref point
LHScoeffderiv = cell(param.nfreqref,param.nthetaref);
RHSderiv = cell(param.nfreqref,param.nthetaref);


for n=1:param.nfreqref
    for m=1:param.nthetaref
        %%%%% LHS %%%%%%
        LHScoeffderiv{n,m}=coeff_derivgen_fun(param.freqref(n),param.thetaref(m));
        %Fix coeff_deriv by replacing all NaN values by 0
        tmpidxnan = find(isnan(LHScoeffderiv{n,m}));
        LHScoeffderiv{n,m}(tmpidxnan) = 0;
        tmpidxinf = find(isinf(LHScoeffderiv{n,m}));
        LHScoeffderiv{n,m}(tmpidxinf) = 0;
        %%%%% RHS %%%%%%
        RHS_tmp = zeros(FEmatrices.size_system,param.nvecfreq,param.nvectheta);
        for ii=1:param.nvecfreq
            for jj=1:param.nvectheta
                %RHS_tmp(:,ii,jj) = build_RHS(param.freqref(n),param.thetaref(m),FEmatrices,[ii,jj],param);
                RHS_tmp(:,ii,jj) = fe_asm(FEmatrices,param,param.freqref(n),param.thetaref(m),[ii,jj]);
            end
        end
        RHSderiv{n,m} = RHS_tmp;
    end
end

end






