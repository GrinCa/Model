function run_MQP()
    fid = fopen('converge.txt','rt');
    while true
        line = fgets(fid);
        if line == -1
            break;
        else
            line = str2num(strtrim(line));
            plot(param.freq,line(2:end),'DisplayName',['ndof = ' num2str(line(1))]);
            hold on
        end
    end
    
    xlabel('Frequency (Hz)');
    ylabel('Mean quadratic pressure (dB)');
    legend();
    hold off
    fclose(fid);
end