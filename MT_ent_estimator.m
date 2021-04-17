function sig = MT_ent_estimator(traj_x_coarse,traj_t_coarse,k_fold,MT_traj)
s_per_pixel = 4.0;

trans_tot = zeros(3);
cond_trans_tot = zeros(3,1);
t_tot = 0;

for i = k_fold
    X = MT_traj(i).traj;
    X = ceil(X(1:traj_t_coarse:end)/traj_x_coarse);
    [trans,cond_trans] = gen_stats(X,1,3);
    trans_tot = trans_tot + trans;
    cond_trans_tot = cond_trans_tot + cond_trans;
    t_tot = t_tot + s_per_pixel*length(X)*traj_t_coarse;
end
t_tot = t_tot/60; % convert to minutes
cond_trans_tot = cond_trans_tot /t_tot; % convert to rate
trans_tot = trans_tot/t_tot; % convert to rate

sig = 0;
 for kk = 1:3
    s = circshift(1:3,kk);
    [~,v] = EntropyEst(trans_tot(s(2),s(3)), trans_tot(s(1),s(2)), trans_tot(s(2),s(1)), cond_trans_tot(s(2)));
    sig = sig + 0.5*v;
   
 end
end