function TL = get_STL(FEmatrices, param, SOL, arg)

[P_inc, P_rad] = get_powers(FEmatrices, param, SOL, arg);

tau_theta = abs(P_rad./P_inc);

tau = zeros(param.nfreq,1);

for ii=1:param.nfreq
    for kk=1:(param.ntheta-1)
        if tau_theta(kk) > 3
            tau_theta(kk) = tau_theta(kk-1);
        end
        tau(ii) = tau(ii) + (tau_theta(kk+1)+tau_theta(kk))*(param.theta(kk+1)+param.theta(kk))/2;
    end    
end

TL = 10*log(1./tau);

end
