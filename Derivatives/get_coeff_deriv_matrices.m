function [coeff_deriv_fun,RHScoeffderiv] = get_coeff_deriv_matrices(derivative_orders,nmat)

deriv_orders_pre_calc = [1,1];
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

coeff_deriv_fun{1,2,1} = @(f,theta) 0;

coeff_deriv_fun{1,2,2} = @(f,theta) 0;

coeff_deriv_fun{2,1,1} = @(f,theta) -4*f^2*pi^2;

coeff_deriv_fun{2,1,2} = @(f,theta) 0;

coeff_deriv_fun{2,2,1} = @(f,theta) -8*f*pi^2;

coeff_deriv_fun{2,2,2} = @(f,theta) 0;

RHScoeffderiv{1,1,1} = @(f,theta,x1,x2) 1000.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*10i)./1711);

RHScoeffderiv{1,1,2} = @(f,theta,x1,x2) -(f.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*10i)./1711).*(x2.*cos(theta) + x1.*sin(theta)).*10000i)./1711;

RHScoeffderiv{1,2,1} = @(f,theta,x1,x2) (pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*10i)./1711).*(x1.*cos(theta) - x2.*sin(theta)).*10000i)./1711;

RHScoeffderiv{1,2,2} = @(f,theta,x1,x2) -(10000.*pi.*exp((f.*pi.*(x1.*cos(theta) - x2.*sin(theta)).*10i)./1711).*(x2.*cos(theta) + x1.*sin(theta)).*(10.*f.*x2.*pi.*sin(theta) - 10.*f.*x1.*pi.*cos(theta) + 1711i))./2927521;

end