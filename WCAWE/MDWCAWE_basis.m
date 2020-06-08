function Wtrans = MDWCAWE_basis(listLHS,coeff_deriv,RHSderiv,nvecfreq,nvectheta)



t_0 = cputime;

ndof = size(listLHS{1},1);
nmatglob = size(listLHS,2);

%warning off all


% nvecfreq = 5;
% nvecsig = 5;

Vtilde = cell(nvecfreq,nvectheta);
V = cell(nvecfreq,nvectheta);
Vmatrix = cell(1,nvectheta);
U = cell(1,nvectheta);



%--------------------------------------------------------------------------
%Initialisation of random cell matrices
%--------------------------------------------------------------------------

% ndof = 216;
% RHSderiv = cell(nvecfreq,nvecsig);
% Kglob = rand(ndof);
% Mglob = rand(ndof);
% listLHS = {Kglob,Mglob};
% %Aglob = rand(ndof,ndof);
% nvecfreq = 5;
% nvecsig = 5;
% nmatglob = 2;
% coeff_deriv = rand(2,nvecfreq,nvecsig);

for jj=1:nvectheta
    U{jj}=[];
    for ii=1:nvecfreq
        Vtilde{ii,jj} = zeros(ndof,1);
        V{ii,jj} = zeros(ndof,1);
%         RHSderiv{ii,jj} = rand(ndof,1);
%         coeff_deriv(1,ii,jj) = 1;
%         coeff_deriv(2,ii,jj) = 1;
    end
end


Aglob = sparse(ndof,ndof);
for kk=1:nmatglob
    Aglob = Aglob + coeff_deriv(kk,1,1)*listLHS{kk};
end

% RHSderiv{1,1} = ones(ndof,1);


%--------------------------------------------------------------------------
% MDWCAWE Algorithm
%--------------------------------------------------------------------------


for jj=1:nvectheta
   if norm(RHSderiv(:,1,jj)) ~= 0
       Vtilde{1,jj} = Aglob\RHSderiv(:,1,jj);
       if jj~=1
           [Vtilde, U] = orthogonalise(1,jj,U,Vtilde,V,nvecfreq);
       end
       U{jj}(1,1) = norm(Vtilde{1,jj});
       V{1,jj} = Vtilde{1,jj}/U{jj}(1,1);
       Vmatrix = concatenate(Vmatrix,V,1,jj);

       for ii=2:nvecfreq
           sumtmp1 = 0;
           sumtmp2 = 0;
           for p1=1:(ii-1)
               PU1 = Pu1(ii,jj,p1,U);
               sumtmp1 = sumtmp1 + RHSderiv(:,p1+1,jj)*PU1(1,ii-p1);  
           end % p1
           for p2=2:(ii-1)
               Aglobderivp = sparse(ndof,ndof);
               for kk=1:nmatglob
                   Aglobderivp = Aglobderivp + coeff_deriv(kk,p2+1,1)*listLHS{kk};
               end % k
               PU2 = Pu2(ii,jj,p2,U);
               sumtmp2 = sumtmp2 + Aglobderivp*Vmatrix{jj}(:,1:ii-p2)*PU2(:,ii-p2);
           end % p2

           Aglobderiv1 = coeff_deriv(1,2,1)*listLHS{1} + coeff_deriv(2,2,1)*listLHS{2};
           Vtilde{ii,jj} = Aglob\(sumtmp1 - Aglobderiv1*V{ii-1,jj} - sumtmp2);
           [Vtilde, U] = orthogonalise(ii,jj,U,Vtilde,V,nvecfreq);
           U{jj}(ii,ii) = norm(Vtilde{ii,jj});
           V{ii,jj} = Vtilde{ii,jj}/U{jj}(ii,ii);
           Vmatrix = concatenate(Vmatrix,V,ii,jj);
       end % i
   end
end % j

Wtrans = full(cell2mat(Vmatrix));


t_end_reduc = cputime-t_0;
outputdisplay = sprintf('[MDWCAWE:INFO] CPUtime for building of WCAWE basis (%d vectors): %.4f s',size(Wtrans,2),t_end_reduc);
disp(outputdisplay);



%--------------------------------------------------------------------------
%Verifs
%--------------------------------------------------------------------------

for ii=1:size(Wtrans,2)
    for jj=1:size(Wtrans,2)
        %disp(sprintf("%d,%d",ii,jj));
        %disp(Wtrans(:,ii)'*Wtrans(:,ii));
    end
end


function [Vtilde, U] = orthogonalise(ii,jj,U,Vtilde,V,nvecfreq)
    for gamma=1:jj
        if gamma==jj
            h=ii-1;
        else
            h=nvecfreq;%ii;%h=nvecfreq;
        end
        for ksi=1:h
            U{gamma}(ksi,ii) = V{ksi,gamma}'*Vtilde{ii,jj}; %
            Vtilde{ii,jj} = Vtilde{ii,jj}-U{gamma}(ksi,ii)*V{ksi,gamma};
        end %ksi
    end %gamma
end


function Pu = Pu1(ii,jj,pp,U)
    Pu = eye(ii-pp);
    for t=1:pp
        Pu = Pu/U{jj}(t:ii-pp+t-1,t:ii-pp+t-1);
    end
end

function Pu = Pu2(ii,jj,pp,U)
    Pu = eye(ii-pp);
    for t=2:pp
        Pu = Pu/U{jj}(t:ii-pp+t-1,t:ii-pp+t-1);
    end
end


function [Vmatrix] = concatenate(Vmatrix,V,ii,jj)
    Vmatrix{jj} = [Vmatrix{jj} V{ii,jj}];
end

end





