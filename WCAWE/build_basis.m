function Wtranssvd = build_basis(FEmatrices,param,nvecfreq,nvectheta,pade_points_freq,pade_points_theta,arg)

Wtrans = [];
for ii=1:length(pade_points_freq)
    for jj=1:length(pade_points_theta)
        if strcmp(arg.algo,'MDWCAWE')
            Wtranstmp = MDWCAWE_basis(FEmatrices.LHS,...
                                      FEmatrices.LHScoeffderiv{pade_points_freq(ii),pade_points_theta(jj)},...
                                      FEmatrices.RHSderiv{pade_points_freq(ii),pade_points_theta(jj)},...
                                      nvecfreq,nvectheta);
            Wtrans = [Wtrans Wtranstmp];
        end
        if strcmp(arg.algo,'WCAWE')
            Wtranstmp_order = param.nvecfreq+param.nvectheta-1;
            LHScoeffderivtmp = zeros(length(FEmatrices.LHS),Wtranstmp_order);
            RHSderivtmp = zeros(FEmatrices.size_system,Wtranstmp_order);
            counter = 1;
            for k = 1:nvecfreq
                if k<nvecfreq
                    LHScoeffderivtmp(:,counter) = FEmatrices.LHScoeffderiv{pade_points_freq(ii),pade_points_theta(jj)}(:,k,1);
                    RHSderivtmp(:,counter) = FEmatrices.RHSderiv{pade_points_freq(ii),pade_points_theta(jj)}(:,k,1);
                    counter = counter+1;
                elseif k==nvecfreq
                    for l = 1:nvectheta
                       LHScoeffderivtmp(:,counter) = FEmatrices.LHScoeffderiv{pade_points_freq(ii),pade_points_theta(jj)}(:,k,l);
                       RHSderivtmp(:,counter) = FEmatrices.RHSderiv{pade_points_freq(ii),pade_points_theta(jj)}(:,k,l);
                       counter = counter+1;
                    end
                end
            end
            Wtranstmp = WCAWE_basis(FEmatrices,LHScoeffderivtmp,RHSderivtmp,Wtranstmp_order);
            Wtrans = [Wtrans Wtranstmp];

            LHScoeffderivtmp = zeros(length(FEmatrices.LHS),Wtranstmp_order);
            RHSderivtmp = zeros(FEmatrices.size_system,Wtranstmp_order);
            counter = 1;
            for l = 1:nvectheta
                if l<nvectheta
                    LHScoeffderivtmp(:,counter) = FEmatrices.LHScoeffderiv{pade_points_freq(ii),pade_points_theta(jj)}(:,1,l);
                    RHSderivtmp(:,counter) = FEmatrices.RHSderiv{pade_points_freq(ii),pade_points_theta(jj)}(:,1,l);
                    counter = counter+1;
                elseif l==nvectheta
                    for k = 1:nvecfreq
                       LHScoeffderivtmp(:,counter) = FEmatrices.LHScoeffderiv{pade_points_freq(ii),pade_points_theta(jj)}(:,k,l);
                       RHSderivtmp(:,counter) = FEmatrices.RHSderiv{pade_points_freq(ii),pade_points_theta(jj)}(:,k,l);
                       counter = counter+1;
                    end % ii
                end
            end % jj
            Wtranstmp = WCAWE_basis(FEmatrices,LHScoeffderivtmp,RHSderivtmp,Wtranstmp_order);
            Wtrans = [Wtrans Wtranstmp(:,2:end-1)];
        end
    end
end


[uu,vv,~] = svd(Wtrans,0);
iiselect = find(diag(vv)>vv(1,1)*1e-15);
Wtranssvd = uu(:,iiselect);
nsvd = size(Wtranssvd,2);
output = sprintf("[SVD:Info] Number of selected vector %d/%d",nsvd,size(Wtrans,2));
disp(output);

end