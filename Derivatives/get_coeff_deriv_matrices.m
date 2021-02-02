function [coeff_deriv_fun,RHScoeffderiv] = get_coeff_deriv_matrices()
coeff_deriv_fun{1,1,1} = @(f,theta) 1 + 1i/20;

coeff_deriv_fun{1,1,2} = @(f,theta) 0;

coeff_deriv_fun{1,2,1} = @(f,theta) 0;

coeff_deriv_fun{1,2,2} = @(f,theta) 0;

coeff_deriv_fun{1,3,1} = @(f,theta) 0;

coeff_deriv_fun{1,3,2} = @(f,theta) 0;

coeff_deriv_fun{1,4,1} = @(f,theta) 0;

coeff_deriv_fun{1,4,2} = @(f,theta) 0;

coeff_deriv_fun{1,5,1} = @(f,theta) 0;

coeff_deriv_fun{1,5,2} = @(f,theta) 0;

coeff_deriv_fun{1,6,1} = @(f,theta) 0;

coeff_deriv_fun{1,6,2} = @(f,theta) 0;

coeff_deriv_fun{1,7,1} = @(f,theta) 0;

coeff_deriv_fun{1,7,2} = @(f,theta) 0;

coeff_deriv_fun{1,8,1} = @(f,theta) 0;

coeff_deriv_fun{1,8,2} = @(f,theta) 0;

coeff_deriv_fun{1,9,1} = @(f,theta) 0;

coeff_deriv_fun{1,9,2} = @(f,theta) 0;

coeff_deriv_fun{1,10,1} = @(f,theta) 0;

coeff_deriv_fun{1,10,2} = @(f,theta) 0;

coeff_deriv_fun{1,11,1} = @(f,theta) 0;

coeff_deriv_fun{1,11,2} = @(f,theta) 0;

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

coeff_deriv_fun{2,6,1} = @(f,theta) 0;

coeff_deriv_fun{2,6,2} = @(f,theta) 0;

coeff_deriv_fun{2,7,1} = @(f,theta) 0;

coeff_deriv_fun{2,7,2} = @(f,theta) 0;

coeff_deriv_fun{2,8,1} = @(f,theta) 0;

coeff_deriv_fun{2,8,2} = @(f,theta) 0;

coeff_deriv_fun{2,9,1} = @(f,theta) 0;

coeff_deriv_fun{2,9,2} = @(f,theta) 0;

coeff_deriv_fun{2,10,1} = @(f,theta) 0;

coeff_deriv_fun{2,10,2} = @(f,theta) 0;

coeff_deriv_fun{2,11,1} = @(f,theta) 0;

coeff_deriv_fun{2,11,2} = @(f,theta) 0;

RHScoeffderiv{1,1,1} = @(f,theta,x1,x2) cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170);

RHScoeffderiv{1,1,2} = @(f,theta,x1,x2) -(f.*pi.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x2.*cos(theta) - x1.*sin(theta)))./170;

RHScoeffderiv{1,2,1} = @(f,theta,x1,x2) -(pi.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)))./170;

RHScoeffderiv{1,2,2} = @(f,theta,x1,x2) - (pi.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x2.*cos(theta) - x1.*sin(theta)))./170 - (f.*pi.^2.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).*(x2.*cos(theta) - x1.*sin(theta)))./28900;

RHScoeffderiv{1,3,1} = @(f,theta,x1,x2) -(pi.^2.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^2)./28900;

RHScoeffderiv{1,3,2} = @(f,theta,x1,x2) (f.*pi.^3.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^2.*(x2.*cos(theta) - x1.*sin(theta)))./4913000 - (pi.^2.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).*(x2.*cos(theta) - x1.*sin(theta)))./14450;

RHScoeffderiv{1,4,1} = @(f,theta,x1,x2) (pi.^3.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^3)./4913000;

RHScoeffderiv{1,4,2} = @(f,theta,x1,x2) (3.*pi.^3.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^2.*(x2.*cos(theta) - x1.*sin(theta)))./4913000 + (f.*pi.^4.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^3.*(x2.*cos(theta) - x1.*sin(theta)))./835210000;

RHScoeffderiv{1,5,1} = @(f,theta,x1,x2) (pi.^4.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^4)./835210000;

RHScoeffderiv{1,5,2} = @(f,theta,x1,x2) (pi.^4.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^3.*(x2.*cos(theta) - x1.*sin(theta)))./208802500 - (f.*pi.^5.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^4.*(x2.*cos(theta) - x1.*sin(theta)))./141985700000;

RHScoeffderiv{1,6,1} = @(f,theta,x1,x2) -(pi.^5.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^5)./141985700000;

RHScoeffderiv{1,6,2} = @(f,theta,x1,x2) - (pi.^5.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^4.*(x2.*cos(theta) - x1.*sin(theta)))./28397140000 - (f.*pi.^6.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^5.*(x2.*cos(theta) - x1.*sin(theta)))./24137569000000;

RHScoeffderiv{1,7,1} = @(f,theta,x1,x2) -(pi.^6.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^6)./24137569000000;

RHScoeffderiv{1,7,2} = @(f,theta,x1,x2) (f.*pi.^7.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^6.*(x2.*cos(theta) - x1.*sin(theta)))./4103386730000000 - (3.*pi.^6.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^5.*(x2.*cos(theta) - x1.*sin(theta)))./12068784500000;

RHScoeffderiv{1,8,1} = @(f,theta,x1,x2) (pi.^7.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^7)./4103386730000000;

RHScoeffderiv{1,8,2} = @(f,theta,x1,x2) (7.*pi.^7.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^6.*(x2.*cos(theta) - x1.*sin(theta)))./4103386730000000 + (f.*pi.^8.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^7.*(x2.*cos(theta) - x1.*sin(theta)))./697575744100000000;

RHScoeffderiv{1,9,1} = @(f,theta,x1,x2) (pi.^8.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^8)./697575744100000000;

RHScoeffderiv{1,9,2} = @(f,theta,x1,x2) (pi.^8.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^7.*(x2.*cos(theta) - x1.*sin(theta)))./87196968012500000 - (f.*pi.^9.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^8.*(x2.*cos(theta) - x1.*sin(theta)))./118587876497000000000;

RHScoeffderiv{1,10,1} = @(f,theta,x1,x2) -(pi.^9.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^9)./118587876497000000000;

RHScoeffderiv{1,10,2} = @(f,theta,x1,x2) - (9.*pi.^9.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^8.*(x2.*cos(theta) - x1.*sin(theta)))./118587876497000000000 - (f.*pi.^10.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^9.*(x2.*cos(theta) - x1.*sin(theta)))./20159939004490000000000;

RHScoeffderiv{1,11,1} = @(f,theta,x1,x2) -(pi.^10.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^10)./20159939004490000000000;

RHScoeffderiv{1,11,2} = @(f,theta,x1,x2) (f.*pi.^11.*sin((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^10.*(x2.*cos(theta) - x1.*sin(theta)))./3427189630763300000000000 - (pi.^10.*cos((f.*pi.*(x1.*cos(theta) + x2.*sin(theta)))./170).*(x1.*cos(theta) + x2.*sin(theta)).^9.*(x2.*cos(theta) - x1.*sin(theta)))./2015993900449000000000;

end