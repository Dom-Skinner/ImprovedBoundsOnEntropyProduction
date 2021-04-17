function [xgs,v] = EntropyEst(nUW, nVU, nUV, nVUPUWV)
% Calculate a lower bound for entropy production given that we have seen
% nVW transitions V -> W, nVU transitions V -> U, nUV transitions U -> V,
% and nVUPUWV transitions W -> V -> U.
% We do this by optimizing over a system with 4 hidden states within the V
% macrostate.

% Sanity checks
assert(isscalar(nUW + nVU + nUV + nVUPUWV),'Error: Non-scalar input')
assert( (nUW >= 0) && (nVU >= 0) && (nUV >= 0) && (nVUPUWV >= 0))

if nVUPUWV == 0
    v = 0;
    xgs = ones(3,1);
    return 
end
min_size = 1e-10; % can't have any input too small for numerical purposes
nUW     = max(nUW,min_size); 
nVU     = max(nVU,min_size);
nUV     = max(nUV,min_size);
nVUPUWV = max(nVUPUWV,min_size);
if min([nUW,nVU,nUV,nVUPUWV]) == min_size % warn user if inputs changed
    warning('Some inputs too small, replaced by default tolerance')
end

n_int = 4; % number of internal states
ObjectiveFunction = @(x) sys_entropy(x,n_int);
nvars = 3*n_int;    % Number of variables
LB = zeros(nvars,1);   % Lower bound
UB = [nUV*ones(n_int,1); nVU*ones(n_int,1); nUW*ones(n_int,1)];  % Upper bound

% Set up the constraints for the problem. Constraints come from probability
% consv. and specifics about the entropy
A = zeros(2 + n_int,3*n_int);
A(1,1:n_int) = 1;
A(2,1+2*n_int:3*n_int) = 1;
for i = 1:n_int
    A(2+i,i) = -1;
    A(2+i,i+n_int) = 1;
    A(2+i,i+2*n_int) = -1;
end
b = [nUV;nUW;zeros(n_int,1)];
Aeq = [ones(1,n_int),-1*ones(1,n_int),zeros(1,n_int)];
beq = nUV-nVU;
ConstraintFunction = @(x) simple_constraint(x,nVUPUWV,n_int); % the one non-linear constraint

% Set up the problem with an initial guess and solve
x0 = [nUV*ones(n_int,1)/n_int;nVU*ones(n_int,1)/n_int;nUW*ones(n_int,1)/n_int];
ConTol = 1e-6;
opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',900,...
                    'MaxFunctionEvaluations',280*nvars,...
                    'FunctionTolerance',1e-7,'ConstraintTolerance',ConTol);
%opts.Display = 'iter';
problem = createOptimProblem('fmincon','objective',ObjectiveFunction,...
            'x0',x0,'lb',LB,'ub',UB,'nonlcon',ConstraintFunction,...
            'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'options',opts);
gs = GlobalSearch('BasinRadiusFactor',0.01);
[xgs,v,EXITFLAG] = run(gs,problem);


% Sometimes optimization fails to converge, due to numerical issues. In
% this case, we rerun with a larger constraint tolerance. Running with a
% large constraint tolerance will only find solutions with a lower
% objective value than the constrained minimum, which is what is needed for
% a lower bound. We warn the user in any case. If this comes up a lot check
% to see if your input parameters are reasonable. A common reason why this
% happens is when you have seen 0 reverse steps vs >1000 forward steps, so
% clearly the statistics are not well sampled.
while EXITFLAG < 1
    ConTol = 4*ConTol;
    opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',800,...
                    'MaxFunctionEvaluations',180*nvars,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
                    'FunctionTolerance',1e-5,'ConstraintTolerance',ConTol);
    problem = createOptimProblem('fmincon','objective',ObjectiveFunction,...
                    'x0',x0,'lb',LB,'ub',UB,'nonlcon',ConstraintFunction,...
                    'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'options',opts);
    [xgs,v,EXITFLAG] = run(gs,problem);
    if (ConTol > 1e-4) && (EXITFLAG < 0)
        warning('Large constraint tolerance had to be used')
    end
end


function [sig,dsig] = sys_entropy(x,n_int)
    % the objective function, in this case the entropy for the system
    niV = x(1:n_int);
    nVi = x(1+n_int:2*n_int);
    niW = x(1+2*n_int:3*n_int);
    nWi = niW + niV - nVi;
    tol_ = 1e-10;
    f = @(x,y) (x - y) .* real(log((x+tol_)./(y+tol_)));
    g = @(x,y) (x - y)./(x+tol_) + real(log((x+tol_)./(y+tol_)));
     
    sig = sum(f(niV , nVi)) + sum( f(niW, nWi));
    
    dsig_dniV = g(niV,nVi);
    dsig_dnVi = g(nVi,niV);
    dsig_dniW = g(niW,nWi);
    dsig = [dsig_dniV;dsig_dnVi;dsig_dniW];
end

function [c, ceq,DC,DCeq] = simple_constraint(x,nVUPUWV,n_int)
    % The one non-linear equality constraint
    niV = x(1:n_int);
    nVi = x(1+n_int:2*n_int);
    niW = x(1+2*n_int:3*n_int);
    tol_ = 1e-10;
    c = [];
    ceq = -sum(nVi.*niW ./(niV+niW + tol_)) + nVUPUWV;
    DC = [];
    dc_dnVi = - niW./(niV + niW + tol_);
    dc_dniV =  (nVi.*niW) ./ (niV + niW + tol_).^2;
    dc_dniW =  -(nVi.*niV) ./ (niV + niW + tol_).^2;
    DCeq = [dc_dniV;dc_dnVi;dc_dniW];
end
end