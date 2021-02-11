function [coeff_deriv_fun,RHScoeffderiv] = get_coeff_deriv_matrices()
coeff_deriv_fun{1,1,1} = @(f,theta) 1;

coeff_deriv_fun{1,1,2} = @(f,theta) 0;

coeff_deriv_fun{1,2,1} = @(f,theta) 0;

coeff_deriv_fun{1,2,2} = @(f,theta) 0;

coeff_deriv_fun{1,3,1} = @(f,theta) 0;

coeff_deriv_fun{1,3,2} = @(f,theta) 0;

coeff_deriv_fun{1,4,1} = @(f,theta) 0;

coeff_deriv_fun{1,4,2} = @(f,theta) 0;

coeff_deriv_fun{1,5,1} = @(f,theta) 0;

coeff_deriv_fun{1,5,2} = @(f,theta) 0;

coeff_deriv_fun{2,1,1} = @(f,theta) -4*f^2*pi^2;

coeff_deriv_fun{2,1,2} = @(f,theta) 0;

coeff_deriv_fun{2,2,1} = @(f,theta) -8*f*pi^2;

coeff_deriv_fun{2,2,2} = @(f,theta) 0;

coeff_deriv_fun{2,3,1} = @(f,theta) -8*pi^2;

coeff_deriv_fun{2,3,2} = @(f,theta) 0;

coeff_deriv_fun{2,4,1} = @(f,theta) 0;

coeff_deriv_fun{2,4,2} = @(f,theta) 0;

coeff_deriv_fun{2,5,1} = @(f,theta) 0;

coeff_deriv_fun{2,5,2} = @(f,theta) 0;

RHScoeffderiv{1,1,1} = @(f,theta,x1,x2) exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170);

RHScoeffderiv{1,1,2} = @(f,theta,x1,x2) (f.*pi.*exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) - x1.*sin(theta)).*1i)./170;

RHScoeffderiv{1,2,1} = @(f,theta,x1,x2) (pi.*exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170;

RHScoeffderiv{1,2,2} = @(f,theta,x1,x2) -(pi.*exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) - x1.*sin(theta)).*(f.*x1.*pi.*cos(theta) + f.*x2.*pi.*sin(theta) - 170i))./28900;

RHScoeffderiv{1,3,1} = @(f,theta,x1,x2) -(pi.^2.*exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) + x2.*sin(theta)).^2)./28900;

RHScoeffderiv{1,3,2} = @(f,theta,x1,x2) -(pi.^2.*exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170).*(f.*x1.*pi.*cos(theta).*1i + f.*x2.*pi.*sin(theta).*1i + 340).*((x2.^2.*sin(2.*theta))./2 - (x1.^2.*sin(2.*theta))./2 + x1.*x2.*cos(2.*theta)))./4913000;

RHScoeffderiv{1,4,1} = @(f,theta,x1,x2) -(pi.^3.*exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) + x2.*sin(theta)).^3.*1i)./4913000;

RHScoeffderiv{1,4,2} = @(f,theta,x1,x2) (pi.^3.*exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) + x2.*sin(theta)).^2.*(x2.*cos(theta) - x1.*sin(theta)).*(f.*x1.*pi.*cos(theta) + f.*x2.*pi.*sin(theta) - 510i))./835210000;

RHScoeffderiv{1,5,1} = @(f,theta,x1,x2) (pi.^4.*exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) + x2.*sin(theta)).^4)./835210000;

RHScoeffderiv{1,5,2} = @(f,theta,x1,x2) (pi.^4.*exp((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) + x2.*sin(theta)).^3.*(x2.*cos(theta) - x1.*sin(theta)).*(f.*x1.*pi.*cos(theta).*1i + f.*x2.*pi.*sin(theta).*1i + 680))./141985700000;

end