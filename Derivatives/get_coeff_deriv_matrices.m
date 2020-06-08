function [coeff_deriv_fun,RHScoeffderiv] = get_coeff_deriv_matrices(derivative_orders,nmat)

deriv_orders_pre_calc = [8,4];
nmat_pre_calc=2;

if nmat > nmat_pre_calc
	disp('[Error] : data mismatch')
	return;
end

for ii=1:length(derivative_orders)
	if derivative_orders(ii)>deriv_orders_pre_calc(ii)
		disp('[Error] : data mismatch')
		return;
	end
end

coeff_deriv_fun{1,1,1} = @(f,theta) 1;

coeff_deriv_fun{1,1,2} = @(f,theta) 0;

coeff_deriv_fun{1,1,3} = @(f,theta) 0;

coeff_deriv_fun{1,1,4} = @(f,theta) 0;

coeff_deriv_fun{1,1,5} = @(f,theta) 0;

coeff_deriv_fun{1,2,1} = @(f,theta) 0;

coeff_deriv_fun{1,2,2} = @(f,theta) 0;

coeff_deriv_fun{1,2,3} = @(f,theta) 0;

coeff_deriv_fun{1,2,4} = @(f,theta) 0;

coeff_deriv_fun{1,2,5} = @(f,theta) 0;

coeff_deriv_fun{1,3,1} = @(f,theta) 0;

coeff_deriv_fun{1,3,2} = @(f,theta) 0;

coeff_deriv_fun{1,3,3} = @(f,theta) 0;

coeff_deriv_fun{1,3,4} = @(f,theta) 0;

coeff_deriv_fun{1,3,5} = @(f,theta) 0;

coeff_deriv_fun{1,4,1} = @(f,theta) 0;

coeff_deriv_fun{1,4,2} = @(f,theta) 0;

coeff_deriv_fun{1,4,3} = @(f,theta) 0;

coeff_deriv_fun{1,4,4} = @(f,theta) 0;

coeff_deriv_fun{1,4,5} = @(f,theta) 0;

coeff_deriv_fun{1,5,1} = @(f,theta) 0;

coeff_deriv_fun{1,5,2} = @(f,theta) 0;

coeff_deriv_fun{1,5,3} = @(f,theta) 0;

coeff_deriv_fun{1,5,4} = @(f,theta) 0;

coeff_deriv_fun{1,5,5} = @(f,theta) 0;

coeff_deriv_fun{1,6,1} = @(f,theta) 0;

coeff_deriv_fun{1,6,2} = @(f,theta) 0;

coeff_deriv_fun{1,6,3} = @(f,theta) 0;

coeff_deriv_fun{1,6,4} = @(f,theta) 0;

coeff_deriv_fun{1,6,5} = @(f,theta) 0;

coeff_deriv_fun{1,7,1} = @(f,theta) 0;

coeff_deriv_fun{1,7,2} = @(f,theta) 0;

coeff_deriv_fun{1,7,3} = @(f,theta) 0;

coeff_deriv_fun{1,7,4} = @(f,theta) 0;

coeff_deriv_fun{1,7,5} = @(f,theta) 0;

coeff_deriv_fun{1,8,1} = @(f,theta) 0;

coeff_deriv_fun{1,8,2} = @(f,theta) 0;

coeff_deriv_fun{1,8,3} = @(f,theta) 0;

coeff_deriv_fun{1,8,4} = @(f,theta) 0;

coeff_deriv_fun{1,8,5} = @(f,theta) 0;

coeff_deriv_fun{1,9,1} = @(f,theta) 0;

coeff_deriv_fun{1,9,2} = @(f,theta) 0;

coeff_deriv_fun{1,9,3} = @(f,theta) 0;

coeff_deriv_fun{1,9,4} = @(f,theta) 0;

coeff_deriv_fun{1,9,5} = @(f,theta) 0;

coeff_deriv_fun{2,1,1} = @(f,theta) -4*f^2*pi^2;

coeff_deriv_fun{2,1,2} = @(f,theta) 0;

coeff_deriv_fun{2,1,3} = @(f,theta) 0;

coeff_deriv_fun{2,1,4} = @(f,theta) 0;

coeff_deriv_fun{2,1,5} = @(f,theta) 0;

coeff_deriv_fun{2,2,1} = @(f,theta) -8*f*pi^2;

coeff_deriv_fun{2,2,2} = @(f,theta) 0;

coeff_deriv_fun{2,2,3} = @(f,theta) 0;

coeff_deriv_fun{2,2,4} = @(f,theta) 0;

coeff_deriv_fun{2,2,5} = @(f,theta) 0;

coeff_deriv_fun{2,3,1} = @(f,theta) -8*pi^2;

coeff_deriv_fun{2,3,2} = @(f,theta) 0;

coeff_deriv_fun{2,3,3} = @(f,theta) 0;

coeff_deriv_fun{2,3,4} = @(f,theta) 0;

coeff_deriv_fun{2,3,5} = @(f,theta) 0;

coeff_deriv_fun{2,4,1} = @(f,theta) 0;

coeff_deriv_fun{2,4,2} = @(f,theta) 0;

coeff_deriv_fun{2,4,3} = @(f,theta) 0;

coeff_deriv_fun{2,4,4} = @(f,theta) 0;

coeff_deriv_fun{2,4,5} = @(f,theta) 0;

coeff_deriv_fun{2,5,1} = @(f,theta) 0;

coeff_deriv_fun{2,5,2} = @(f,theta) 0;

coeff_deriv_fun{2,5,3} = @(f,theta) 0;

coeff_deriv_fun{2,5,4} = @(f,theta) 0;

coeff_deriv_fun{2,5,5} = @(f,theta) 0;

coeff_deriv_fun{2,6,1} = @(f,theta) 0;

coeff_deriv_fun{2,6,2} = @(f,theta) 0;

coeff_deriv_fun{2,6,3} = @(f,theta) 0;

coeff_deriv_fun{2,6,4} = @(f,theta) 0;

coeff_deriv_fun{2,6,5} = @(f,theta) 0;

coeff_deriv_fun{2,7,1} = @(f,theta) 0;

coeff_deriv_fun{2,7,2} = @(f,theta) 0;

coeff_deriv_fun{2,7,3} = @(f,theta) 0;

coeff_deriv_fun{2,7,4} = @(f,theta) 0;

coeff_deriv_fun{2,7,5} = @(f,theta) 0;

coeff_deriv_fun{2,8,1} = @(f,theta) 0;

coeff_deriv_fun{2,8,2} = @(f,theta) 0;

coeff_deriv_fun{2,8,3} = @(f,theta) 0;

coeff_deriv_fun{2,8,4} = @(f,theta) 0;

coeff_deriv_fun{2,8,5} = @(f,theta) 0;

coeff_deriv_fun{2,9,1} = @(f,theta) 0;

coeff_deriv_fun{2,9,2} = @(f,theta) 0;

coeff_deriv_fun{2,9,3} = @(f,theta) 0;

coeff_deriv_fun{2,9,4} = @(f,theta) 0;

coeff_deriv_fun{2,9,5} = @(f,theta) 0;

RHScoeffderiv{1,1,1} = @(f,theta,x1,x2) exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170);

RHScoeffderiv{1,1,2} = @(f,theta,x1,x2) -(f.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*1i)./170;

RHScoeffderiv{1,1,3} = @(f,theta,x1,x2) - (f.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170 - (f.^2.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./28900;

RHScoeffderiv{1,1,4} = @(f,theta,x1,x2) (f.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(510.*f.*x2.*pi.*sin(theta) - 510.*f.*x1.*pi.*cos(theta) + f.^2.*x2.^2.*pi.^2.*cos(theta).^2.*1i + f.^2.*x1.^2.*pi.^2.*sin(theta).^2.*1i + f.^2.*x1.*x2.*pi.^2.*cos(theta).*sin(theta).*2i + 28900i))./4913000;

RHScoeffderiv{1,1,5} = @(f,theta,x1,x2) (f.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170 + (f.^2.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./7225 - (3.*f.^2.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./28900 + (f.^4.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4)./835210000 + (f.^3.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*3i)./2456500;

RHScoeffderiv{1,2,1} = @(f,theta,x1,x2) (pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170;

RHScoeffderiv{1,2,2} = @(f,theta,x1,x2) -(pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(f.*x2.*pi.*sin(theta) - f.*x1.*pi.*cos(theta) + 170i))./28900;

RHScoeffderiv{1,2,3} = @(f,theta,x1,x2) (f.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./28900 - (f.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./14450 - (pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170 - (f.^2.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./4913000;

RHScoeffderiv{1,2,4} = @(f,theta,x1,x2) (pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(202300.*f.*x2.*pi.*sin(theta) - 202300.*f.*x1.*pi.*cos(theta) - (f.^3.*x1.^3.*pi.^3.*cos(theta))./4 + (f.^3.*x2.^3.*pi.^3.*sin(theta))./4 - f.^2.*x1.^2.*pi.^2.*cos(2.*theta).*510i + f.^2.*x2.^2.*pi.^2.*cos(2.*theta).*510i + (f.^3.*x1.^3.*pi.^3.*cos(3.*theta))./4 + (f.^3.*x2.^3.*pi.^3.*sin(3.*theta))./4 + f.^2.*x1.*x2.*pi.^2.*sin(2.*theta).*1020i + (f.^3.*x1.^2.*x2.*pi.^3.*sin(theta))./4 - (3.*f.^3.*x1.*x2.^2.*pi.^3.*cos(3.*theta))./4 - (3.*f.^3.*x1.^2.*x2.*pi.^3.*sin(3.*theta))./4 - (f.^3.*x1.*x2.^2.*pi.^3.*cos(theta))./4 + 4913000i))./835210000;

RHScoeffderiv{1,2,5} = @(f,theta,x1,x2) (pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170 - (f.^2.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*3i)./4913000 + (f.^3.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4)./208802500 + (2.*f.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./7225 - (7.*f.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./28900 + (f.^2.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*11i)./2456500 + (f.^4.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./141985700000 - (3.*f.^3.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./417605000;

RHScoeffderiv{1,3,1} = @(f,theta,x1,x2) -(pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./28900;

RHScoeffderiv{1,3,2} = @(f,theta,x1,x2) (pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(f.*x1.*pi.*cos(theta).*1i - f.*x2.*pi.*sin(theta).*1i + 340).*((x1.^2.*sin(2.*theta))./2 - (x2.^2.*sin(2.*theta))./2 + x1.*x2.*cos(2.*theta)))./4913000;

RHScoeffderiv{1,3,3} = @(f,theta,x1,x2) (pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./14450 - (pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./14450 + (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./4913000 + (f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./835210000 - (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./1228250;

RHScoeffderiv{1,3,4} = @(f,theta,x1,x2) (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*3i)./2456500 - (2.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)))./7225 - (3.*f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)))./417605000 + (3.*f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^3)./835210000 - (f.^3.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^2.*1i)./141985700000 - (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^2.*13i)./4913000;

RHScoeffderiv{1,3,5} = @(f,theta,x1,x2) (2.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./7225 - (2.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./7225 + (3.*f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4)./208802500 + (3.*f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^4)./835210000 - (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*13i)./4913000 + (f.^3.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./17748212500 - (f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./20880250 - (f.^3.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*3i)./70992850000 - (f.^4.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^2)./24137569000000 + (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*13i)./1228250;

RHScoeffderiv{1,4,1} = @(f,theta,x1,x2) -(pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./4913000;

RHScoeffderiv{1,4,2} = @(f,theta,x1,x2) (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^2.*(f.*x2.*pi.*sin(theta) - f.*x1.*pi.*cos(theta) + 510i))./835210000;

RHScoeffderiv{1,4,3} = @(f,theta,x1,x2) (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*3i)./4913000 - (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*3i)./2456500 - (f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^4)./835210000 + (3.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./417605000 + (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./141985700000;

RHScoeffderiv{1,4,4} = @(f,theta,x1,x2) (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*3i)./2456500 - (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^2.*21i)./4913000 + (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^4.*3i)./141985700000 - (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^2.*9i)./141985700000 + (f.^3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^3)./24137569000000 - (9.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)))./417605000 + (19.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^3)./835210000;

RHScoeffderiv{1,4,5} = @(f,theta,x1,x2) (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^5.*3i)./141985700000 - (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*21i)./4913000 + (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*3i)./245650 + (3.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4)./104401250 + (19.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^4)./835210000 - (33.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./208802500 + (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).*9i)./35496425000 - (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*29i)./70992850000 + (3.*f.^3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^4)./12068784500000 - (3.*f.^3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^2)./6034392250000 - (f.^4.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./4103386730000000;

RHScoeffderiv{1,5,1} = @(f,theta,x1,x2) (pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^4)./835210000;

RHScoeffderiv{1,5,2} = @(f,theta,x1,x2) -(pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^3.*(f.*x1.*pi.*cos(theta).*1i - f.*x2.*pi.*sin(theta).*1i + 680))./141985700000;

RHScoeffderiv{1,5,3} = @(f,theta,x1,x2) (3.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./208802500 - (pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^4)./208802500 - (f.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^5.*1i)./141985700000 + (f.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./17748212500 - (f.^2.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^4)./24137569000000;

RHScoeffderiv{1,5,4} = @(f,theta,x1,x2) (pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^3)./20880250 - (3.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)))./104401250 - (f.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^2.*9i)./35496425000 - (3.*f.^2.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^5)./24137569000000 + (3.*f.^2.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^3)./6034392250000 + (f.^3.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^4.*1i)./4103386730000000 + (f.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^4.*1i)./5679428000;

RHScoeffderiv{1,5,5} = @(f,theta,x1,x2) (3.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4)./104401250 + (pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^4)./20880250 - (3.*f.^2.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^6)./24137569000000 - (12.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./52200625 + (f.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^5.*1i)./5679428000 - (f.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*31i)./17748212500 + (19.*f.^2.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^4)./6034392250000 - (9.*f.^2.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^2)./3017196125000 + (f.^3.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^5.*3i)./2051693365000000 - (f.^3.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./256461670625000 + (f.^4.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^4)./697575744100000000 + (f.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).*3i)./4437053125;

RHScoeffderiv{1,6,1} = @(f,theta,x1,x2) (pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^5.*1i)./141985700000;

RHScoeffderiv{1,6,2} = @(f,theta,x1,x2) -(pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^4.*(f.*x2.*pi.*sin(theta) - f.*x1.*pi.*cos(theta) + 850i))./24137569000000;

RHScoeffderiv{1,6,3} = @(f,theta,x1,x2) (pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./7099285000 - (pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^5.*1i)./28397140000 + (f.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^6)./24137569000000 - (f.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^4)./2413756900000 - (f.^2.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^5.*1i)./4103386730000000;

RHScoeffderiv{1,6,4} = @(f,theta,x1,x2) (pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^4.*13i)./28397140000 - (pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^2.*3i)./7099285000 + (3.*f.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^3)./1206878450000 - (f.^2.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^6.*3i)./4103386730000000 + (f.^2.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^4.*3i)./820677346000000 - (f.^3.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^5)./697575744100000000 - (31.*f.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^5)./24137569000000;

RHScoeffderiv{1,6,5} = @(f,theta,x1,x2) (pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^5.*13i)./28397140000 - (f.^2.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^7.*3i)./4103386730000000 + (pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).*3i)./3549642500 - (pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*11i)./3549642500 - (31.*f.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^6)./24137569000000 + (f.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^4)./60343922500 - (3.*f.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^2)./301719612500 + (f.^2.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^5.*47i)./2051693365000000 - (f.^2.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^3.*3i)./102584668250000 - (3.*f.^3.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^6)./348787872050000000 + (f.^3.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^4)./34878787205000000 + (f.^4.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^5.*1i)./118587876497000000000;

RHScoeffderiv{1,7,1} = @(f,theta,x1,x2) -(pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^6)./24137569000000;

RHScoeffderiv{1,7,2} = @(f,theta,x1,x2) (pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^5.*(f.*x1.*pi.*cos(theta).*1i - f.*x2.*pi.*sin(theta).*1i + 1020))./4103386730000000;

RHScoeffderiv{1,7,3} = @(f,theta,x1,x2) (3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^6)./12068784500000 - (3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^4)./2413756900000 + (f.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^7.*1i)./4103386730000000 - (f.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^5.*3i)./1025846682500000 + (f.^2.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^6)./697575744100000000;

RHScoeffderiv{1,7,4} = @(f,theta,x1,x2) (3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^3)./603439225000 - (3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^5)./754299031250 + (f.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^4.*9i)./410338673000000 + (3.*f.^2.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^7)./697575744100000000 - (9.*f.^2.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^5)./348787872050000000 - (f.^3.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^6.*1i)./118587876497000000000 - (f.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^6.*37i)./4103386730000000;

RHScoeffderiv{1,7,5} = @(f,theta,x1,x2) (3.*f.^2.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^8)./697575744100000000 - (3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^6)./754299031250 + (21.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^4)./603439225000 - (9.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^2)./603439225000 - (f.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^7.*37i)./4103386730000000 + (f.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^5.*147i)./1025846682500000 - (f.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^3.*3i)./25646167062500 - (7.*f.^2.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^6)./43598484006250000 + (9.*f.^2.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^4)./34878787205000000 - (f.^3.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^7.*3i)./59293938248500000000 + (f.^3.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^5.*3i)./14823484562125000000 - (f.^4.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^6)./20159939004490000000000;

RHScoeffderiv{1,8,1} = @(f,theta,x1,x2) -(pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^7.*1i)./4103386730000000;

RHScoeffderiv{1,8,2} = @(f,theta,x1,x2) (pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^6.*(f.*x2.*pi.*sin(theta) - f.*x1.*pi.*cos(theta) + 1190i))./697575744100000000;

RHScoeffderiv{1,8,3} = @(f,theta,x1,x2) (pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^7.*7i)./4103386730000000 - (pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^5.*21i)./2051693365000000 - (f.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^8)./697575744100000000 + (7.*f.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^6)./348787872050000000 + (f.^2.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^7.*1i)./118587876497000000000;

RHScoeffderiv{1,8,4} = @(f,theta,x1,x2) (pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^4.*21i)./410338673000000 - (pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^6.*133i)./4103386730000000 - (63.*f.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^5)./348787872050000000 + (f.^2.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^8.*3i)./118587876497000000000 - (f.^2.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^6.*21i)./118587876497000000000 + (f.^3.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^7)./20159939004490000000000 + (43.*f.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^7)./697575744100000000;

RHScoeffderiv{1,8,5} = @(f,theta,x1,x2) (f.^2.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^9.*3i)./118587876497000000000 - (pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^7.*133i)./4103386730000000 + (pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^5.*21i)./60343922500000 - (pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^3.*21i)./102584668250000 + (43.*f.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^8)./697575744100000000 - (203.*f.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^6)./174393936025000000 + (21.*f.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^4)./17439393602500000 - (f.^2.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^7.*13i)./11858787649700000000 + (f.^2.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^5.*63i)./29646969124250000000 + (3.*f.^3.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^8)./10079969502245000000000 - (7.*f.^3.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^6)./5039984751122500000000 - (f.^4.*pi.^11.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^7.*1i)./3427189630763300000000000;

RHScoeffderiv{1,9,1} = @(f,theta,x1,x2) (pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^8)./697575744100000000;

RHScoeffderiv{1,9,2} = @(f,theta,x1,x2) -(pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^7.*(f.*x1.*pi.*cos(theta).*1i - f.*x2.*pi.*sin(theta).*1i + 1360))./118587876497000000000;

RHScoeffderiv{1,9,3} = @(f,theta,x1,x2) (7.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^6)./87196968012500000 - (pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^8)./87196968012500000 - (f.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^9.*1i)./118587876497000000000 + (f.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^7.*1i)./7411742281062500000 - (f.^2.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^8)./20159939004490000000000;

RHScoeffderiv{1,9,4} = @(f,theta,x1,x2) (11.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^7)./43598484006250000 - (21.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^5)./43598484006250000 - (f.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^6.*21i)./14823484562125000000 - (3.*f.^2.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^9)./20159939004490000000000 + (3.*f.^2.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^7)./2519992375561250000000 + (f.^3.*pi.^11.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^8.*1i)./3427189630763300000000000 + (f.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^8.*49i)./118587876497000000000;

RHScoeffderiv{1,9,5} = @(f,theta,x1,x2) (11.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^8)./43598484006250000 - (3.*f.^2.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^10)./20159939004490000000000 - (7.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^6)./2179924200312500 + (21.*pi.^8.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^4)./8719696801250000 + (f.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^9.*49i)./118587876497000000000 - (f.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^7.*67i)./7411742281062500000 + (f.*pi.^9.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^5.*21i)./1852935570265625000 + (37.*f.^2.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^8)./5039984751122500000000 - (21.*f.^2.*pi.^10.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^6)./1259996187780625000000 + (f.^3.*pi.^11.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^9.*3i)./1713594815381650000000000 - (f.^3.*pi.^11.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^7.*1i)./107099675961353125000000 + (f.^4.*pi.^12.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^8)./582622237229761000000000000;

end