function show_graph(arg,data,mesh,param)

if strcmp(arg.type,'compare_results')
    sizemesh = load('sizemesh.txt');
    if length(sizemesh) > 1
        disp('There are several sizemesh');
        disp('Please choose one amoung the following :');
        for iii=1:length(sizemesh)
            disp([num2str(iii),' -> sizemesh = ',num2str(sizemesh(iii))]);
        end
    end
    id_sample =  ['_sizemesh_',num2str(sizemesh)];
    SOLFE =      struct2cell(load(['Matrices/',mesh.file,'/',param.path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat']));
    SOLMDWCAWE = struct2cell(load(['Matrices/',mesh.file,'/',param.path2,'/SOLMDWCAWE',id_sample,'.mat']));
    SOLWCAWE =   struct2cell(load(['Matrices/',mesh.file,'/',param.path2,'/SOLWCAWE',id_sample,'.mat']));
    
    SOLFE = SOLFE{1};
    SOLMDWCAWE = SOLMDWCAWE{1};
    SOLWCAWE = SOLWCAWE{1};
    
    
    
    if ~isempty(find(size(SOLFE)==1)) || length(size(SOLFE)) == 2
        Norm_SOLFE = zeros(param.nfreq,param.ntheta);
        Norm_SOLMDWCAWE = zeros(param.nfreq,param.ntheta);
        Norm_SOLWCAWE = zeros(param.nfreq,param.ntheta);
        if length(size(SOLFE)) == 2
            for ii=1:size(SOLFE,2)
                Norm_SOLFE(ii) = norm(SOLFE(:,ii));
                Norm_SOLMDWCAWE(ii) = norm(SOLMDWCAWE(:,ii));
                Norm_SOLWCAWE(ii) = norm(SOLWCAWE(:,ii));
            end
        else
            for ii=1:size(SOLFE,3)
                Norm_SOLFE(ii) = norm(SOLFE(:,1,ii));
                Norm_SOLMDWCAWE(ii) = norm(SOLMDWCAWE(:,1,ii));
                Norm_SOLWCAWE(ii) = norm(SOLWCAWE(:,1,ii));
            end
        end
        show2d(SOLFE,SOLMDWCAWE,SOLWCAWE,param);
    else
        Norm_SOLFE = zeros(param.nfreq,param.ntheta);
        Norm_SOLMDWCAWE = zeros(param.nfreq,param.ntheta);
        Norm_SOLWCAWE = zeros(param.nfreq,param.ntheta);
        for ii=1:size(SOLFE,2)
            for jj=1:size(SOLFE,3)
                Norm_SOLFE(ii,jj) = norm(SOLFE(:,ii,jj));
                Norm_SOLMDWCAWE(ii,jj) = norm(SOLMDWCAWE(:,ii,jj));
                Norm_SOLWCAWE(ii,jj) = norm(SOLWCAWE(:,ii,jj));
            end
        end
        show3d(Norm_SOLFE,Norm_SOLMDWCAWE,Norm_SOLWCAWE,param,mesh);
    end

elseif strcat(arg.type,'plotTL')
    
    if ~isempty(find(size(data{1})==1))
        show2d(data,param,mesh,arg)
    else
        show3d(data,param,mesh,arg)
    end
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
        ylabel("TL");
        title("");
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
            title_label = arg.title;
            title(title_label);
        end
    else
        figure
        for ii=1:length(VALUES)
            surf(X,Y,VALUES{ii}','DisplayName',arg.label{ii});
            xlabel('freq');
            ylabel('theta');
            zlabel(arg.zlabel);
            title_label = arg.title;
            title(title_label);
            hold on
        end
    end
    
    argtmp.gcf = gcf;
    argtmp.type = 'save_plot';
    argtmp.save_name = [arg.name_plot '_' param.interval_detail_str '.png'];
    IO_data(argtmp,param,mesh);

end


