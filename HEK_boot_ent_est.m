function [sig_lb,sig_ub,sig_med] = HEK_boot_ent_est(traj_arr_in,ub,lb,t_tot_post,NBoot)
k = size(traj_arr_in,2);
sig_boot = zeros(NBoot,1);
for i = 1:length(sig_boot)
    r = randi([1 k],1,k);
    sig_boot(i) = HEK_ent_est(traj_arr_in,r,ub,lb,t_tot_post);
end
sig_lb  = quantile(sig_boot,0.025);
sig_ub  = quantile(sig_boot,0.975);
sig_med = quantile(sig_boot,0.5);
end