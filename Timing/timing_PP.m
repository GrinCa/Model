function timing_PP(timing, param, arg)

if strcmp(arg, 'save')
    save_timing(timing, param)
end

if strcmp(arg, 'display')
    display(param)
end

end

function save_timing(timing, param)
keys_cell = timing.keys();
values_cell = timing.values();
try
    data = struct2cell(load(['Matrices/',param.filename,'/timing.mat']));
    timing_backup = data{1};
catch
    timing_backup = containers.Map;
end

for ii=1:length(keys_cell)
    timing_backup(keys_cell{ii}) = values_cell{ii};
end

save(['Matrices/', param.filename, '/timing.mat'], 'timing_backup');

end


function display(param)

try
    data = struct2cell(load(['Matrices/',param.filename,'/timing.mat']));
    timing = data{1};
catch
    disp("Can't display times");
end

keys_cell = timing.keys();
values_cell = timing.values();

for ii=1:length(keys_cell)
    disp([keys_cell{ii}, '  ---->  ',num2str(values_cell{ii}), ' secs']);
end

end
