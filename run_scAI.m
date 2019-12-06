function result = run_scAI(X1,X2,Ks,alpha,lambda,gamma,s,Inits,repeat,stop_rule,seeds)
if ~exist('alpha','var') || isempty(alpha)
    alpha = 1;
end

if ~exist('lambda','var') || isempty(lambda)
    lambda = 10000;
end

if ~exist('gamma','var') || isempty(gamma)
    gamma = 1;
end

if ~exist('s','var') || isempty(s)
    s = 0.25;
end

if ~exist('stop_rule','var') || isempty(stop_rule)
    stop_rule = 1;
end

if ~exist('repeat','var') || isempty(repeat)
    repeat = 10;
end

if ~exist('Inits','var')
    Inits = [];
end

if ~exist('seeds','var')
    seeds = 1:repeat;
end

n = length(Ks); result = cell(n,repeat);
for i = 1:n
    K = Ks(i);
    for j = 1:repeat
        if ~isempty(Inits)
            [W1,W2,H,Z,R] = scAI(X1,X2,K,alpha,lambda,gamma,s,stop_rule,seeds(j),'W1_INIT',Inits{i,j}.W1,'W2_INIT',Inits{i,j}.W2,'H_INIT',Inits{i,j}.H,'Z_INIT',Inits{i,j}.Z,'R',Inits{i,j}.R);
        else
            [W1,W2,H,Z,R] = scAI(X1,X2,K,alpha,lambda,gamma,s,stop_rule,seeds(j));
        end
        result{i,j}.W1 = W1; result{i,j}.W2 = W2; result{i,j}.H = H; result{i,j}.Z = Z; result{i,j}.R = R;
    end
end