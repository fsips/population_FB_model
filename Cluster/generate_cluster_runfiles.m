function [] = generate_cluster_runfiles()

%% Control plots
f_calls = fopen(['Cluster_calls.m'], 'w');

fprintf(f_calls, ['\n CONTROL CALLS \n']);

n = 10;
for it = [1 2 3 4 5 6 7 8 11 12 13 14 15 16 17 18 19 20 21 22 23 24]; 
    for fxr = 0:1
        % ADD QUOTES
        call = generate_cluster_call(['setDirs; generate_control_plot(',num2str(fxr),', ',num2str(it),', ',num2str(n),')'], ['LOG_CONTROL_', num2str(fxr),'_',num2str(it),'_',num2str(n)]);
        fprintf(f_calls, [call, ' \n']);
    end
end
fclose(f_calls);

