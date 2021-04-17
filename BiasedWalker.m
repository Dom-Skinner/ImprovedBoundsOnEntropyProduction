% Example of a biased random walker that can switch directions. If we only
% see the trajectory of the walker, then every forward path is as likely as
% every reverse path, leading to most estimators to conclude no entropy is
% being produced. Here we use our estimator to get a non-trivial bound on
% the entropy production rate.

% r the probability of changing bias directions, p,q the bias. Here, we
% solve for a range of p to illustrate the problem.
r = 0.05;
p = linspace(r/2,1-r-0.05,20);
q = 1 - r - p;

% Estimate the entropy production rate for the biased walker for each
% (p,q,r) value
sig = zeros(size(p));
for i = 1:length(p)
    
    N = 800000;
    [X,t] = generate_walker_traj(p(i),q(i),p(i),q(i),r,N); %generate a trajectory

    [trans,cond_trans] = gen_stats(X,t(end),3); % calculate the observables

    [~,v1] = EntropyEst(trans(3,1), trans(1,2), trans(2,1), cond_trans(2));
    [~,v2] = EntropyEst(trans(1,2), trans(2,3), trans(3,2), cond_trans(3));
    [~,v3] = EntropyEst(trans(2,3), trans(3,1), trans(1,3), cond_trans(1));
    sig(i) = 0.5*(v1+v2+v3);
    
end

% Now calculate the exact solution
p_fine = linspace(r/2,1-r-0.05,100);
q_fine = 1 - r - p_fine;
sig_true = (p_fine-q_fine).*log(p_fine./q_fine);


plot(p_fine,sig_true)
hold on
scatter(p,sig)
legend('Exact solution','Estimator')
xlabel('p')
ylabel('Rate of entropy production')


function [X,t] = generate_walker_traj(p1,q1,p2,q2,r,N)
    % Generates a trajectory of length N of the biased random walker.
    % To create a continuous time markov chain, simulate a discrete time
    % MC and then draw waiting times from an exponential distribution
    P = [0  q1 p1 r  0  0; ...
         p1 0  q1 0  r  0; ...
         q1 p1 0  0  0  r; ...
         r  0  0  0  p2 q2; ...
         0  r  0  q2 0  p2; ...
         0  0  r  p2 q2 0];
    mc = dtmc(P);
    X = simulate(mc,N);
    dt  = exprnd(1,N,1);
    t = [0; cumsum(dt)];
end

function OW = observed_walk(X)
    % Converts the MC X on 6 states into a walk on Z 
    n = 3;
    Y = mod(X,n); 
    OW = zeros(size(X));
    OW(1) = X(1);
    for i = 2:length(Y)
        if (Y(i) == 0) && (Y(i-1) == n-1)
            OW(i) = OW(i-1) + 1;
        elseif (Y(i) == n-1) && (Y(i-1) == 0)
            OW(i) = OW(i-1) - 1;
        else
            OW(i) = OW(i-1) + Y(i) - Y(i-1);
        end
    end
end