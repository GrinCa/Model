function timing = run_STL(param, mesh, timing, sizemesh)

    clear FEmatrices SOLFE SOLWCAWE SOLMDWCAWE;
    arg.sizemesh = sizemesh;
    arg.type = 'preload';
    FEmatrices = IO_data(arg,param,mesh);
    
    VALUES_SOL = cell(1,3); % index 2
    VALUES_name = cell(1,3);
    % is for MDWCAWE and WCAWE
    exist_FE = 1;
    exist_MDWCAWE = 1;
    exist_WCAWE = 1;
    error_loading_FE = 0;
    error_loading_MDWCAWE = 0;
    error_loading_WCAWE = 0;
    
    %----------------------------------------------------------------------
    % Compute FE solution
    %----------------------------------------------------------------------
    
    try
        path = ['STL/',mesh.file,'/Sizemesh[',num2str(sizemesh),']/',param.path1,'/FE_SOL.mat'];
        if exist(path)
            tmp = load(path);
            VALUES_SOL{1} = tmp.FE_STL;
            VALUES_name{1} = 'FE';
            exist_FE = 0;
        end
        if exist_FE
            arg.type = 'load_FE';
            SOLFE = IO_data(arg,param,mesh);
        end
    catch
        disp("SOLFE not find for the given parameters");
        error_loading_FE = 1;
        exist_FE = 0;
    end
    if exist_FE
        arg.type = 'calculateTL';
        start = cputime;
        VALUES_SOL{1} = post_process(FEmatrices,param,arg,SOLFE);
        timing(['FE_', param.study, '_', param.path1, '_STL']) = cputime - start;
        VALUES_name{1} = 'FE';
        FE_STL = VALUES_SOL{1};
        save(['STL/',mesh.file,'/Sizemesh[',num2str(sizemesh),']/',param.path1,'/FE_SOL.mat'],'FE_STL');
    end
    if error_loading_FE
        VALUES_SOL{1} = zeros(param.nfreq,param.ntheta);
        VALUES_name{1} = 'FE';
    end
    
    %----------------------------------------------------------------------
    % Compute MDWCAWE solution
    %----------------------------------------------------------------------

    try
        path = ['STL/',mesh.file,'/Sizemesh[',num2str(sizemesh),']/',param.path2,'/FE_MDWCAWE.mat'];
        if exist(path)
            tmp = load(path);
            VALUES_SOL{2} = tmp.MDWCAWE_STL;
            VALUES_name{2} = 'MDWCAWE';
            exist_MDWCAWE = 0;
        end
        if exist_MDWCAWE
            arg.type = 'load_MDWCAWE';
            SOLMDWCAWE = IO_data(arg,param,mesh);
        end
    catch
        disp("SOLMDWCAWE not find for the given parameters");
        error_loading_MDWCAWE = 1;
        exist_MDWCAWE = 0;
    end
    if exist_MDWCAWE
        arg.type = 'calculateTL';
        start = cputime;
        VALUES_SOL{2} = post_process(FEmatrices,param,arg,SOLMDWCAWE{1});
        VALUES_name{2} = 'MDWCAWE';
        timing(['MDWCAWE_', param.study, '_', param.path2, '_STL']) = cputime - start;
        MDWCAWE_STL = VALUES_SOL{2};
        save(['STL/',mesh.file,'/Sizemesh[',num2str(sizemesh),']/',param.path2,'/FE_MDWCAWE.mat'],'MDWCAWE_STL');
    end
    if error_loading_MDWCAWE
        VALUES_SOL{2} = zeros(param.nfreq,param.ntheta);
        VALUES_name{2} = 'MDWCAWE';
    end
    
    %----------------------------------------------------------------------
    % Compute MDWCAWE solution
    %----------------------------------------------------------------------    
    
    try
        path = ['STL/',mesh.file,'/Sizemesh[',num2str(sizemesh),']/',param.path2,'/FE_WCAWE.mat'];
        if exist(path)
            tmp = load(path);
            VALUES_SOL{3} = tmp.WCAWE_STL;
            VALUES_name{3} = 'WCAWE';
            exist_WCAWE = 0;
        end
        if exist_WCAWE
            arg.type = 'load_WCAWE';
            SOLWCAWE = IO_data(arg,param,mesh);
        end
    catch
        disp("SOLWCAWE not find for the given parameters");
        error_loading_MDWCAWE = 1;
        exist_WCAWE = 0;
    end
    
    if exist_WCAWE
        start = cputime;
        arg.type = 'calculateTL';
        VALUES_SOL{3} = post_process(FEmatrices,param,arg,SOLWCAWE{1});
        VALUES_name{3} = 'WCAWE';
        timing(['WCAWE_', param.study, '_', param.path2, '_STL']) = cputime - start;
        WCAWE_STL = VALUES_SOL{3};
        save(['STL/',mesh.file,'/Sizemesh[',num2str(sizemesh),']/',param.path2,'/FE_WCAWE.mat'],'WCAWE_STL');
    end
    if error_loading_MDWCAWE
        VALUES_SOL{3} = zeros(param.nfreq,param.ntheta);
        VALUES_name{3} = 'WCAWE';
    end
   
    argcomp.ylabel = 'TL (dB)';
    argcomp.title = 'Comparison FE MDWCAWE';
    argcomp.type = 'plotTL';
    argcomp.split = 0;
    argcomp.name_plot = 'Comparison_FE_MDWCAWE';
    argcomp.label = VALUES_name(1:2);
    argcomp.external_plot.is_needed = false;
    argcomp.path = param.path2;
    show_graph(argcomp,VALUES_SOL(1:2),mesh,param);
    
    argcomp.ylabel = 'TL (dB)';
    argcomp.title = 'Comparison FE WCAWE';
    argcomp.type = 'plotTL';
    argcomp.split = 0;
    argcomp.name_plot = 'Comparison_FE_WCAWE';
    argcomp.label = VALUES_name(1:2:end);
    argcomp.external_plot.is_needed = false;
    argcomp.path = param.path2;
    show_graph(argcomp,VALUES_SOL(1:2:end),mesh,param);
    
    argtmp.type = 'rel_error';
    argtmp.REF_SOLUTION = VALUES_SOL(1);
    argtmp.APPROX_SOLUTION = VALUES_SOL(2:end);
    argcomp.type = 'plotTL';
    argcomp.ylabel = 'Relative Error (Log_{10})';
    argcomp.title = '';
    argcomp.split = 1;
    argcomp.name_plot = 'Relative_error';
    argcomp.label = VALUES_name(2:end);
    argcomp.external_plot.is_needed = false;
    argcomp.path = param.path2;
    show_graph(argcomp,post_process(FEmatrices,param,argtmp),mesh,param);

end