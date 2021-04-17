function [sig1,sig2] = entropy_sweep(Path,n_bins,T)
logfn = @(x,y) (x-y)*log(x/y); % for the naive estimator

sig1 = zeros(size(n_bins));
sig2 = zeros(size(n_bins));

for i = 1:length(n_bins)
     X =  floor(3*n_bins(i)*Path);
    [trans,cond_trans] = gen_stats(X,T,3);
    for kk = 1:3
        s = circshift(1:3,kk);
        ntrials = 4;
        z = zeros(ntrials,1);
        for k = 1:ntrials
        [~,v] = EntropyEst(trans(s(2),s(3)), trans(s(1),s(2)), trans(s(2),s(1)), cond_trans(s(2)));
        z(k) = v;
        end
        sig2(i) = sig2(i) + 0.5*min(z);
        sig1(i) = sig1(i) + logfn(trans(s(1),s(2)),trans(s(2),s(1))); 
    end
end
end