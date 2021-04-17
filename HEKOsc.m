% For figure 5. Adjusting the lb and ub parameters slightly adjusts the
% results. They were not optimized over to avoid overfitting this small
% dataset. Also, more consistent bootstrap results come with more samples 
% (NBoot~200 >>1), but this will make the code slow. 
clear
load('HEK_data')

lb = 0.39; 
ub = 0.57;
lb_200 = 0.52;
ub_200 = 0.65;
NBoot = 50;

[sig_lb_pre,sig_ub_pre,sig_med_pre] = HEK_boot_ent_est(traj_arr_pre,ub,lb,t_tot_pre,NBoot);
[sig_lb_100,sig_ub_100,sig_med_100] = HEK_boot_ent_est(traj_arr_post_100,ub, lb,t_tot_post,NBoot);
[sig_lb_200,sig_ub_200,sig_med_200] = HEK_boot_ent_est(traj_arr_post_200,ub_200,lb_200,t_tot_post,NBoot);

sig_med = 60*[sig_med_pre;sig_med_100;sig_med_200];
sig_lb = 60*[sig_lb_pre;sig_lb_100;sig_lb_200];
sig_ub = 60*[sig_ub_pre;sig_ub_100;sig_ub_200];
errorbar(1:3,sig_med,sig_med-sig_lb,sig_ub-sig_med,...
    'o','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
ylim([0,15])
