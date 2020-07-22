function SOLMDWCAWE = Solve_MDWCAWE(FEmatrices,param,MDWCAWE,interval_freq,interval_theta)


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

SOLMDWCAWE = zeros(length(param.idx_out),length(interval_freq),length(interval_theta));

for ii=1:length(interval_freq)
   for jj=1:length(interval_theta)
      Aglob_red = sparse(size(LHS_red{1},1),size(LHS_red{1},2));
      for kk=1:nmat
         Aglob_red = Aglob_red + FEmatrices.LHScoeffderiv_fun{kk,1,1}(interval_freq(ii),interval_theta(jj))*LHS_red{kk};
      end %kk
      %%% Calculation of RHS @ (f,theta) %%%
      RHS = build_RHS(interval_freq(ii),interval_theta(jj),FEmatrices,[1,1],param);
      RHS = sparse(MDWCAWE'*RHS(:));
      solp = Aglob_red\RHS;
      resp_P = MDWCAWE*solp;
      SOLMDWCAWE(:,ii,jj) = resp_P(param.idx_out);
   end % jj
end % ii


end

