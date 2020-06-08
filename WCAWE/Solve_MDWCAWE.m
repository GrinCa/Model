function SOLMDWCAWE = Solve_MDWCAWE(FEmatrices,MDWCAWE,param)


nmat = length(FEmatrices.LHS);
LHS_red = cell(1,length(FEmatrices.LHS));
%--------------------------------------------------------------------------
% Project matrices on WCAWE basis and setup RHS
%--------------------------------------------------------------------------
for ii = 1:nmat
   LHS_red{ii} = sparse(MDWCAWE'* FEmatrices.LHS{ii}*MDWCAWE);
end %ii

%--------------------------------------------------------------------------
% Frequency/Fi loops calculation
%--------------------------------------------------------------------------

SOLMDWCAWE = zeros(length(param.idx_out),param.nfreq,param.ntheta);

for ii=1:param.nfreq
   for jj=1:param.ntheta
      Aglob_red = sparse(size(LHS_red{1},1),size(LHS_red{1},2));
      for kk = 1:nmat
         Aglob_red = Aglob_red + FEmatrices.LHScoeffderiv_fun{kk,1,1}(param.freq(ii),param.theta(jj))*LHS_red{kk};
      end %kk
      %%% Calculation of RHS @ (f,theta) %%%
      RHS = build_RHS(param.freq(ii),param.theta(jj),FEmatrices,FEmatrices.RHScoeffderiv_fun,[1,1],param.direction);
      RHS = sparse(MDWCAWE'*RHS(:));
      solp = Aglob_red\RHS;
      resp_P = MDWCAWE*solp;
      SOLMDWCAWE(:,ii,jj) = resp_P(param.idx_out);
   end % jj
end % ii


end

