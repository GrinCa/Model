function normalize_error()

    VALUES = cell(1,2);
    arg.sizemesh = sizemesh;
    arg.type = 'preload';
    FEmatrices = IO_data(arg,param,mesh);
    arg.type = 'load_FE';
    arg.REF_SOLUTION = IO_data(arg,param,mesh);
    arg.type = 'load_MDWCAWE';
    arg.APPROX_SOLUTION = IO_data(arg,param,mesh);
    arg.type = 'normalize_error';
    VALUES{1} = log10(post_process(FEmatrices,param,arg));
    arg.split = false;
    arg.label = {'MDWCAWE','WCAWE'};
    arg.zlabel = 'error';
    arg.title = {'FE/MDWCAWE','FE/WCAWE'};
    arg.external_plot.is_needed = false;
    arg.name_plot = 'relative error norm';
    %show_graph(arg,VALUES,mesh,param);
    
    
    arg.type = 'load_WCAWE';
    arg.APPROX_SOLUTION = IO_data(arg,param,mesh);
    arg.type = 'normalize_error';
    arg.split = true;
    VALUES{2} = log10(post_process(FEmatrices,param,arg));
    show_graph(arg,VALUES,mesh,param);

end