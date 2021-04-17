function sig = HEK_ent_est(traj_arr,boot_arr,ub,lb,t_tot)

trans_tot = zeros(3);
cond_trans_tot = zeros(3,1);

for i = boot_arr

    traj = traj_arr(:,i);
    traj(traj<lb) = 0;
    traj(traj> ub) = 2;
    traj( (traj~=0) & (traj ~=2)) = 1;

    [trans,cond_trans] = gen_stats(traj,1,3);
     trans_tot = trans_tot + trans;
    cond_trans_tot = cond_trans_tot + cond_trans;
end
trans_tot = trans_tot/(t_tot*length(boot_arr));
cond_trans_tot = cond_trans_tot/(t_tot*length(boot_arr));

[~,sig] = EntropyEst(trans_tot(1,2), trans_tot(3,1), trans_tot(1,3), cond_trans_tot(1));
%{
sig = 0;
for kk = 1:3
    s = circshift(1:3,kk);
    [~,v] = EntropyEst(trans_tot(s(2),s(3)), trans_tot(s(2),s(1)), trans_tot(s(1),s(2)), cond_trans_tot(s(2)));
    sig = sig + 0.5*v;
end
 %}
end