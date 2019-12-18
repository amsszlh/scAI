function [W1,W2,H,Z,R] = scAI(X1,X2,K,s,alpha,lambda,gamma,stop_rule,seed,varargin)
% This algorithm is to solve:
% min alpha*||X1-W1*H||_F^2+||X2*(Z.*R)-W2H||_F^2+lambda*||Z-H'H||_F^2+gamma*Sum(||H(:,j)||_1^2)
% s.t. W1, W2, H,Z >= 0
% if lambda = 0, let gamma=0, Z=eye(n),R = ones(n),the objective function is degenerated into joint NMF  
% <Inputs>
%         X1: scRNA-seq data matrix
%         X2: scATAC-seq/single cell DNA methylation data matrix         
%         K: the number of clusters
%         s: usually defined from 0.2 to 0.3, which is used to control aggregation numbers
%         alpha: the parameter to balance information trusted for X1 and X2
%         lambda: the parameter to balance the error of Z and H
%         gamma: the parameter to control sparsity of each row of H
%         stop_rule: If stop_rule = 1 (default), the algorithm stops if the
%         iter high than 500. If stop_value = 2, the algorithm stops if it
%         satisfies (objs(iter) - objs(iter+1))/objs(1)< 10^(-6).
%         seed: The seed number used to generate repeatable result
%         (Below are optional arguments: can be set by providing name-value pairs)
%         R:  rand aggregation among subpopulation R = binornd(1,s,size(X1,2),size(X1,2));
%         W1_INIT : (p x K) initial value for W1.
%         W2_INIT : (q x K) initial value for W2.
%         H_INIT : (K x n) initial value for H.
%         Z_INIT : (n x n) initial value for Z.
% <Outputs>
%         W1, W2: the basis nonnegative low rank matrix
%         H: the coefficient nonnegative low rank matrix
%         Z: the similarity matrix
%         R: rand aggregation matrix
%% Initialization W1,W2,H and Z
[p,n] = size(X1); q = size(X2,1);
rng(seed)
W1 = rand(p,K); W2 = rand(q,K); H = rand(K,n); Z = rand(n); R = binornd(1,s,n,n);
Maxiter = 500; 
%% Read optional parameters
if (rem(length(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for i = 1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'R',            R = varargin{i+1};
            case 'W1_INIT',      W1 = varargin{i+1};
            case 'W2_INIT',      W2 = varargin{i+1};
            case 'H_INIT',       H = varargin{i+1};
            case 'Z_INIT',       Z = varargin{i+1};
            case 'MAXITER',      Maxiter = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end
if isempty(W1)
    W1 = rand(p,K);
end
if isempty(W2)
    W2 = rand(q,K);
end
if isempty(H) 
    H = rand(K,n);
end
if isempty(Z)
    Z = rand(n);
end
if isempty(R)
    R = binornd(1,s,n,n);
end
% whether lambda = 0 or not
if lambda == 0
    gamma = 0;
	R = ones(n); Z = eye(n);
end
X2 = sparse(X2); 
XtX2 = full(X2'*X2); obj_old = 1; Index = R == 0;
for iter = 1:Maxiter
     % normalized H
    lib = sum(H,2); H = H./repmat(lib,1,n);
    % update W1
    HHt = H*H'; X1Ht = X1*H'; W1HHt = W1*HHt;
    W1 = W1.*X1Ht./(W1HHt+eps);
    % update W2
    ZR = Z; ZR(Index) = 0; 
    ZRHt = ZR*H'; X2ZRHt = X2*ZRHt; W2HHt = W2*HHt;
    W2 = W2.*X2ZRHt./(W2HHt+eps);
    % update H
    W1tX1 = W1'*X1; W2tX2 = W2'*X2; W2tX2ZR = W2tX2*ZR; HZZt = H*(Z+Z'); W1tW1 = W1'*W1; W2tW2 = W2'*W2;
    H = H.*(alpha*W1tX1+W2tX2ZR+lambda*HZZt)./((alpha*W1tW1+W2tW2+2*lambda*HHt+gamma*ones(K))*H+eps);
    % update Z
    if lambda ~= 0
        HtH = H'*H;
        RX2tW2H = (W2tX2)'*H; RX2tW2H(Index) = 0;  
        XtX2ZRR = XtX2*ZR; XtX2ZRR(Index) = 0;
        Z = Z.*(RX2tW2H+lambda*HtH)./(XtX2ZRR+lambda*Z+eps);
    end
    if stop_rule == 2
        obj = norm(X1-W1*H,'fro')^2 + norm(X2*(Z.*R)-W2*H,'fro')^2 + lambda*norm(Z-HtH,'fro')^2 + gamma*sum(sum(H).^2);
        if ((obj_old-obj)/obj_old < 10^(-6) && iter > 1) || iter == Maxiter
            break;
        end
        obj_old = obj;
    end
end

