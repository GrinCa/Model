function [coeff_deriv_fun,RHScoeffderiv] = get_coeff_deriv_matrices(derivative_orders,nmat)

deriv_orders_pre_calc = [3,4];
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

RHScoeffderiv{1,1,1} = @(f,theta,x1,x2) exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170)./800;

RHScoeffderiv{1,1,2} = @(f,theta,x1,x2) -(f.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*1i)./136000;

RHScoeffderiv{1,1,3} = @(f,theta,x1,x2) - (f.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./136000 - (f.^2.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./23120000;

RHScoeffderiv{1,1,4} = @(f,theta,x1,x2) (f.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(510.*f.*x2.*pi.*sin(theta) - 510.*f.*x1.*pi.*cos(theta) + f.^2.*x2.^2.*pi.^2.*cos(theta).^2.*1i + f.^2.*x1.^2.*pi.^2.*sin(theta).^2.*1i + f.^2.*x1.*x2.*pi.^2.*cos(theta).*sin(theta).*2i + 28900i))./3930400000;

RHScoeffderiv{1,1,5} = @(f,theta,x1,x2) (f.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./136000 + (f.^2.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./5780000 - (3.*f.^2.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./23120000 + (f.^4.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4)./668168000000 + (f.^3.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*3i)./1965200000;

RHScoeffderiv{1,2,1} = @(f,theta,x1,x2) (pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./136000;

RHScoeffderiv{1,2,2} = @(f,theta,x1,x2) -(pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(f.*x2.*pi.*sin(theta) - f.*x1.*pi.*cos(theta) + 170i))./23120000;

RHScoeffderiv{1,2,3} = @(f,theta,x1,x2) (f.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./23120000 - (f.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./11560000 - (pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./136000 - (f.^2.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./3930400000;

RHScoeffderiv{1,2,4} = @(f,theta,x1,x2) (pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(202300.*f.*x2.*pi.*sin(theta) - 202300.*f.*x1.*pi.*cos(theta) - (f.^3.*x1.^3.*pi.^3.*cos(theta))./4 + (f.^3.*x2.^3.*pi.^3.*sin(theta))./4 - f.^2.*x1.^2.*pi.^2.*cos(2.*theta).*510i + f.^2.*x2.^2.*pi.^2.*cos(2.*theta).*510i + (f.^3.*x1.^3.*pi.^3.*cos(3.*theta))./4 + (f.^3.*x2.^3.*pi.^3.*sin(3.*theta))./4 + f.^2.*x1.*x2.*pi.^2.*sin(2.*theta).*1020i + (f.^3.*x1.^2.*x2.*pi.^3.*sin(theta))./4 - (3.*f.^3.*x1.*x2.^2.*pi.^3.*cos(3.*theta))./4 - (3.*f.^3.*x1.^2.*x2.*pi.^3.*sin(3.*theta))./4 - (f.^3.*x1.*x2.^2.*pi.^3.*cos(theta))./4 + 4913000i))./668168000000;

RHScoeffderiv{1,2,5} = @(f,theta,x1,x2) (pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).*1i)./136000 - (f.^2.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*3i)./3930400000 + (f.^3.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4)./167042000000 + (f.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./2890000 - (7.*f.*pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./23120000 + (f.^2.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*11i)./1965200000 + (f.^4.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./113588560000000 - (3.*f.^3.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./334084000000;

RHScoeffderiv{1,3,1} = @(f,theta,x1,x2) -(pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./23120000;

RHScoeffderiv{1,3,2} = @(f,theta,x1,x2) (pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(f.*x1.*pi.*cos(theta).*1i - f.*x2.*pi.*sin(theta).*1i + 340).*((x1.^2.*sin(2.*theta))./2 - (x2.^2.*sin(2.*theta))./2 + x1.*x2.*cos(2.*theta)))./3930400000;

RHScoeffderiv{1,3,3} = @(f,theta,x1,x2) (pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./11560000 - (pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./11560000 + (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./3930400000 + (f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./668168000000 - (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./982600000;

RHScoeffderiv{1,3,4} = @(f,theta,x1,x2) (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*3i)./1965200000 - (pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)))./2890000 - (3.*f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)))./334084000000 + (3.*f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^3)./668168000000 - (f.^3.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^2.*1i)./113588560000000 - (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^2.*13i)./3930400000;

RHScoeffderiv{1,3,5} = @(f,theta,x1,x2) (pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2)./2890000 - (pi.^2.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^2)./2890000 + (3.*f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4)./167042000000 + (3.*f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^4)./668168000000 - (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*13i)./3930400000 + (f.^3.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./14198570000000 - (f.^2.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./16704200000 - (f.^3.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*3i)./56794280000000 - (f.^4.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^2)./19310055200000000 + (f.*pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*13i)./982600000;

RHScoeffderiv{1,4,1} = @(f,theta,x1,x2) -(pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./3930400000;

RHScoeffderiv{1,4,2} = @(f,theta,x1,x2) (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^2.*(f.*x2.*pi.*sin(theta) - f.*x1.*pi.*cos(theta) + 510i))./668168000000;

RHScoeffderiv{1,4,3} = @(f,theta,x1,x2) (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*3i)./3930400000 - (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*3i)./1965200000 - (f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^4)./668168000000 + (3.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./334084000000 + (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./113588560000000;

RHScoeffderiv{1,4,4} = @(f,theta,x1,x2) (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*3i)./1965200000 - (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^2.*21i)./3930400000 + (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^4.*3i)./113588560000000 - (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^2.*9i)./113588560000000 + (f.^3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)).^3)./19310055200000000 - (9.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^3.*(x1.*cos(theta) - x2.*sin(theta)))./334084000000 + (19.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).*(x1.*cos(theta) - x2.*sin(theta)).^3)./668168000000;

RHScoeffderiv{1,4,5} = @(f,theta,x1,x2) (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^5.*3i)./113588560000000 - (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^3.*21i)./3930400000 + (pi.^3.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).*3i)./196520000 + (3.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4)./83521000000 + (19.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x1.*cos(theta) - x2.*sin(theta)).^4)./668168000000 - (33.*f.*pi.^4.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^2)./167042000000 + (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).*9i)./28397140000000 - (f.^2.*pi.^5.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^3.*29i)./56794280000000 + (3.*f.^3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^2.*(x1.*cos(theta) - x2.*sin(theta)).^4)./9655027600000000 - (3.*f.^3.*pi.^6.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^2)./4827513800000000 - (f.^4.*pi.^7.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*1i)./170).*(x2.*cos(theta) + x1.*sin(theta)).^4.*(x1.*cos(theta) - x2.*sin(theta)).^3.*1i)./3282709384000000000;

end