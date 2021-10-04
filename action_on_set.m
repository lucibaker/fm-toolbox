addpath('Z:\luci\MATLAB\large-particles')
run_id_list = {'run-e' 'run-f' 'run-g' 'run-h'};

for i = 3:length(run_id_list)
    run_id = run_id_list{i};
    disp(run_id)
    cd(['../' run_id])
%     disp('img_to_bmp')
%     img_to_bmp(run_id)
    disp('print_sphere_mask')
    print_sphere_mask(run_id)
end
