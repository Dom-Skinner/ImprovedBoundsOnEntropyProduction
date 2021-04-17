function [trans,cond_trans] = gen_stats(X,T,n)
% Given a Markov chain with trajectory X in total time T, we record the
% rate of transitions A->B, as well as the rate of conditional transitions,
% A -> B -> C.
% In theory this code will work for any discrete X, making it periodic with
% period n. In practice n=3 is what has been used

% make periodic
X = mod(X,n);
X(X==0) = n;

% Since we do not use the time at which things happened, make X so that
% each value is different from the previous one
%{
i = 1;
while i < length(X)
    if X(i) == X(i+1)
        X(i+1) = [];
    else
        i = i+ 1;
    end
end
%}
not_repeat = X(2:end) ~= X(1:end-1);
X = X([true;not_repeat]);

% Count the rate at which we see transitions i -> j
trans = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            trans(i,j) = sum((X(1:end-1) == i) & (X(2:end) == j))/T;
        end
    end
end
    
% Count the rate at which we see transitions i-1 -> i -> i+1
cond_trans = zeros(n,1);
cond_trans(1) = sum((X(1:end-2) == n) & (X(2:end-1) == 1) & (X(3:end) == 2))/T;
cond_trans(n) = sum((X(1:end-2) == n-1) & (X(2:end-1) == n) & (X(3:end) == 1))/T;
for i = 2:n-1
    cond_trans(i) = sum((X(1:end-2) == i-1) & (X(2:end-1) == i) & (X(3:end) == i+1))/T;
end
end
