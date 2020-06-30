function [Wtrans,Ucoeff] = WCAWE_basis(FEmatrices,coeff_deriv,RHSderiv,nbasevec)
% Note that Aglob should be factorized!



t_0 = cputime;
nderiv = size(RHSderiv,2);
nmat = length(FEmatrices.LHS);

% Setup MUMPS for solution alone (multiple RHS)

%--------------------------------------------------------------------------
% Tests on number of basis vector requested, and update size of Wtrans and Ucoeff
%--------------------------------------------------------------------------
if (nderiv<nbasevec)
   disp('Error, too many basis vectors requested!');
   return;
else
end

% nbasevec_current = size(Wtrans,2);

% if (nbasevec_current==0)
%    Ucoeff = zeros(nbasevec,nbasevec);
%    Wtrans = zeros(size(RHSderiv,1),nbasevec);
% 
% elseif (nbasevec_current<nbasevec)
%    Ucoeff(nbasevec_current+1:nbasevec,:) = zeros(nbasevec-nbasevec_current,nbasevec_current);
%    Ucoeff(:,nbasevec_current+1:nbasevec) = zeros(nbasevec,nbasevec-nbasevec_current);
% 
%    Wtrans(:,nbasevec_current+1:nbasevec) = zeros(size(Wtrans,1),nbasevec-nbasevec_current);
% elseif (nbasevec_current==nbasevec) % nothing will be done then...
% else % (nbasevec_current>nbasevec)
%    outputdisplay = sprintf('   WARNING: fewer basis vectors asked than existing!');
%    disp(outputdisplay)
% end



%--------------------------------------------------------------------------
% Calculate or complete WCAWE basis
%--------------------------------------------------------------------------

% calculation of Aglob = K-omega^2*M
Aglob = sparse(FEmatrices.size_system,FEmatrices.size_system);
for kk=1:nmat
    Aglob = Aglob+coeff_deriv(kk,1,1)*FEmatrices.LHS{kk};
end

%initialization of Mumps instance
id=initmumps();
id=zmumps(id);
id.JOB=4;
id=zmumps(id,Aglob);

Wtrans = [];

for kk=1:nbasevec
   % Sum associated with RHS vector
    if kk==1
        RHS = RHSderiv(:,kk);
    else
        RHS = sparse(size(RHSderiv(:,1),1),1);
    end % if
    for jj=1:kk-1
        Pu = eye(kk-jj,kk-jj);
        for ll=1:jj
            Ublock = Ucoeff(ll:kk-jj+ll-1,ll:kk-jj+ll-1);
            Pu = Ublock\Pu;
        end  % ll
        RHS = RHS + RHSderiv(:,jj+1)*Pu(1,kk-jj);
    end % jj
   
    % Sum associated with recursive terms
    for jj=1:kk-1 % Order of derivative for system matrix: (jj-1)
        Pu = eye(kk-jj,kk-jj);
        for ll=2:jj
             Ublock = Ucoeff(ll:kk-jj+ll-1,ll:kk-jj+ll-1);
             Pu = Ublock\Pu;
        end % ll
        for ll=1:nmat
            RHS = RHS - coeff_deriv(ll,jj+1)*FEmatrices.LHS{ll}*Wtrans(:,1:kk-jj)*Pu(:,kk-jj);
        end % ll
    end % jj
   
    % Modified Gram-Schmidt
    id.RHS = RHS;
    id.JOB = 3;
    id = zmumps(id,Aglob);
    Wtrans(:,kk) = id.SOL;
    for jj=1:kk-1
        Ucoeff(jj,kk) = Wtrans(:,jj)'*Wtrans(:,kk);
        Wtrans(:,kk) = Wtrans(:,kk)-Ucoeff(jj,kk)*Wtrans(:,jj);
    end % jj
    Ucoeff(kk,kk) = norm(Wtrans(:,kk));
    Wtrans(:,kk) = Wtrans(:,kk)/Ucoeff(kk,kk);
end % kk

id.JOB = -2;
id = zmumps(id,Aglob);

t_end_reduc = cputime-t_0;
outputdisplay = sprintf('[WCAWE:INFO] CPUtime for building of WCAWE basis (%d vectors): %.4f s',nbasevec,t_end_reduc);
disp(outputdisplay);

