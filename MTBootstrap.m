% For Figure 4
clear
load('MT_traj')

Nboot = 100; % For high quality results this needs to be large, but will take a long time to run.
x_coarse = 1:4;
t_coarse = 2:4;

N = length(x_coarse) + length(t_coarse);
sig_lb = zeros(N,1);
sig_ub = zeros(N,1);
sig_med = zeros(N,1);

for j = 1:N
    sig_boot = zeros(Nboot,1);
    for i = 1:length(sig_boot)
        r = randi([1 81],1,81);
        if j > length(x_coarse)
            sig_boot(i) = MT_ent_estimator(t_coarse(j-length(x_coarse)),1,r,MT_traj);
        else
            sig_boot(i) = MT_ent_estimator(1,x_coarse(j),r,MT_traj);
        end
    end
    sig_lb(j)  = quantile(sig_boot,0.025);
    sig_ub(j)  = quantile(sig_boot,0.975);
    sig_med(j) = quantile(sig_boot,0.5);
end
errorbar(1:4,sig_med(1:4),sig_med(1:4)-sig_lb(1:4),sig_ub(1:4)-sig_med(1:4),...
    'o','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
xlabel('Spatial coarse graining')
ylabel('Entropy production rate (k_B min^-1)')

figure
sig_med_t = [sig_med(1);sig_med(5:end)];
sig_lb_t =  [sig_lb(1);sig_lb(5:end)];
sig_ub_t =  [sig_ub(1);sig_ub(5:end)];
errorbar(1:4,sig_med_t,sig_med_t-sig_lb_t,sig_ub_t-sig_med_t,...
    'o','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
xlabel('Temporal coarse graining')
ylabel('Entropy production rate (k_B min^-1)')
