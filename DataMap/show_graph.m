function show_graph(arg,data,mesh,param)

if ~isempty(find(size(data{1})==1))
    show2d(data,param,mesh,arg)
else
    show3d(data,param,mesh,arg)
end
    
end


function show2d(VALUES,param,mesh,arg)
    X = {param.freq,param.theta};
    sizesubset = [param.nfreq,param.ntheta];
    xlabeltxt = {"freq",'theta'};
    
    figure
    for ii=1:length(VALUES)
        plot(X{find(sizesubset>1)},VALUES{ii},'DisplayName',arg.label{ii});
        xlabel(xlabeltxt{find(sizesubset>1)});
        ylabel(arg.ylabel);
        title(arg.title{ii});
        legend;
        hold on
    end
    
    argtmp.gcf = gcf;
    argtmp.type = 'save_plot';
    argtmp.save_name = [arg.name_plot '_' param.interval_detail_str '.png'];
    IO_data(argtmp,param,mesh);

end


function show3d(VALUES,param,mesh,arg)
    
    [X,Y] = meshgrid(param.freq,param.theta);
    
    if arg.split
        for ii=1:length(VALUES)
            figure
            surf(X,Y,VALUES{ii}','DisplayName',arg.label{ii});
            xlabel('freq');
            ylabel('theta');
            zlabel(arg.zlabel);
            title_label = arg.title{ii};
            title(title_label);
        end
    else
        figure
        for ii=1:length(VALUES)
            surf(X,Y,VALUES{ii}','DisplayName',arg.label{ii});
            xlabel('freq');
            ylabel('theta');
            zlabel(arg.zlabel);
            title_label = arg.title{ii};
            title(title_label);
            hold on
        end
        if arg.external_plot.is_needed
            for ii=1:length(arg.external_plot.VALUES)
                surf(X,Y,arg.external_plot.VALUES{ii}','DisplayName',arg.external_plot.label{ii});
            end
        end
    end
    
    argtmp.gcf = gcf;
    argtmp.type = 'save_plot';
    argtmp.save_name = [arg.name_plot '_' param.interval_detail_str '.png'];
    IO_data(argtmp,param,mesh);

end


