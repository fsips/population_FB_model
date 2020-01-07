function [runstring] = generate_cluster_call(name1, name2)
runstring = ['nohup matlab -nojvm -r ',name1,' -logfile out',name2,'.OUT </dev/null &'];
