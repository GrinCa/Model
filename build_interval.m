function param = build_interval(param)

%build sub_interval with the different freqref
index_freqref = zeros(1,param.nfreqref);
index_thetaref = zeros(1,param.nthetaref);
for ii=1:param.nfreqref
    [~,index_freqref(ii)] = min(abs(param.freq-param.freqref(ii)));
end
for ii=1:param.nthetaref
    [~,index_thetaref(ii)] = min(abs(param.theta-param.thetaref(ii)));
end

index_ref = {index_freqref;
             index_thetaref};

param.sub_interval = {cell(1,length(param.interval_construct{1}));  %freq
                      cell(1,length(param.interval_construct{2}))}; %theta

id_cut = { zeros(1,length(param.interval_construct{1})-1);
           zeros(1,length(param.interval_construct{2})-1) }; %freq + angle

for n=1:2
    for ii=1:length(param.interval_construct{n})-1
        if param.interval_construct{n}{ii}(end) == param.interval_construct{n}{ii+1}(1)
            id_cut{n}(ii) = index_ref{n}(param.interval_construct{n}{ii}(end));
        else
            id_cut{n}(ii) = int16(( index_ref{n}(param.interval_construct{n}{ii}(end)) + ...
                                    index_ref{n}(param.interval_construct{n}{ii+1}(1)) )/2);
        end
    end
end



caracteristic_index = {[0 id_cut{1} param.nfreq];
                       [0 id_cut{2} param.ntheta]};
                   
interval = {param.freq;
            param.theta};
        
increment = {param.freqincr;
             param.thetaincr};
         
% interval_index is needed for the reconstruction of the solution on the
% whole interval
param.interval_index = { cell(1,length(param.interval_construct{1}));  %freq
                         cell(1,length(param.interval_construct{2})) }; %theta
         
for n=1:2
    if isempty(id_cut(n,:))
        param.sub_interval{n} = param.freq;
    else
        for ii=1:(length(caracteristic_index{n})-1)
            left_edge = caracteristic_index{n}(ii)+1;
            right_edge = caracteristic_index{n}(ii+1);
            param.interval_index{n}{ii} = left_edge : 1 : right_edge;
            param.sub_interval{n}{ii} = interval{n}(left_edge) : increment{n} : interval{n}(right_edge);
        end
    end
end

end

