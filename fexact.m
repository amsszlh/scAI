function pval = fexact( varargin )
%FEXACT Fisher's exact test
%
% Fisher's exact test is a statistical test used to compare binary outcomes
% between two groups. For example, a laboratory test may be devised to
% predict outcomes in clinical trials. The question is whether the outcomes
% are more favorable (e.g. the patient Responds to treatment) in the
% subjects testing postive than in the group overall. In the 2x2 table
% below the test results are shown in the rows and the number of patients
% responding (R), or not-responding are shown in columns.  
%        R       NR
%  pos   a        c       K
%  neg   b        d       -
% total  N        -       M
%
% a,b,c,d are the numbers of positive and negative results for the cases
% and controls. K, N, M are the row, column and grand totals.
%
% This function returns a pvalue against the null hypothesis that a
% responders out of K sampled is likely if the population sample has N out
% of M. A small p-value is taken as evidence that a/K differs from N/M.
%
%
% usage
%   p = fexact(X,y)
%       y is a vetor of status (1=case/0=control). X is a MxP matrix of
%       binary results (pos=1/neg=0). P can be very large for
%       genotyping assays, each of which can be considered a different way
%       to categorize the cases and controls.
%
%  p = fexact( a,M,K,N, options)
%      This version is intended for use when you start with a contingency
%      matrix rather than raw test results and outcomes. N is the number of
%      positive outcomes (scalar). M is the total number of observations. a
%      and K are P-vectors. No checks for valid input are made.
%
%
%  p = fexact(..., 'options', values)
%      options: 
%      test  'f', 'l'. Sets the test to be either f(exact) or
%            l(iebermeister). Chris Rorden added this functionality for
%            single-sided tests. 
%      tail  l(eft)','r(ight)', 'b(oth)'
%            left tail tests a negative association. p( x<=a, M, K, N)
%            right tail tests a positive association. p( x<=b, M, K, M-N)
%            both (default) is either a positive or negative association.
%            this is the sum from the left and right tails of the
%            distribution including all x where p(x|M,K,N) is less than
%            or equal p(a|M,K,N).
%      perm  Q, where Q is an integer. Permutation testing when X has multiple Columns. 
%            to correct for multiple testing. This is necessary if some of
%            the test results are correlated (e.g. whole genome genotype
%            assays). The entire set of tests is repeated Q times
%            with permuted y variables. The reported p-values are corrected
%            for multiple tests by interpolating into the emprical
%            distribution of the minimum p-value obtained from each Q*repsz
%            rounds . To use this option the X and y calling convention
%            must be used (as opposed to (a,M,K,N,...) convention. and Q
%            must be greater than 1.
%      repsz S, where S is an integer. Used with the PERM option to specify 
%            the block size if Q is large (default=100). 
%            Larger S is faster but requires more memory. 
%
%
% NOTES and LIMITATIONS:
%      This function is extremely fast when doing multiple tests and is
%      acceptable with a small number of tests. However, it uses large
%      tracks of memory. On my 2Gb home Intel Core2 Quad Core Vista machine
%      this function does 250,000 tests with 100 observations each in 0.10
%      seconds. I run out of memory when I do more
%
% example
%  x = unidrnd(2, 200, 1000)-1;         % simulate test results
%  y = unidrnd(2,200,1)-1;              % simulate known response rates
%  p = fexact( x, y );                  % generates p-values for 1000 tests in x
%  p1 = fexact( x, y, 'perm', 10); % reports p-value relative to empirial
%                                               % cdf, since these are randomly generated 
%                                               % the values in p1 should be near  1
%   
% References:
%   Seneta E, Phipps MC (2001) On the Comparison of Two Observed Frequencies. Biometrical Journal, 43(1):23?43, 2001.
%   Phipps MC (2003)Inequalities between hypergeometric tails. Journal of Applied Mathematics and Decision Sciences. 7(3): 165-174
%   http://www.emis.de/journals/HOA/ADS/Volume7_3/174.pdf
%
% $Id: fexact.m 990 2014-12-12 16:38:19Z mboedigh $
% Copyright 2012 Mike Boedigheimer
% Amgen Inc.
% Department of Computational Biology
% mboedigh@amgen.com

p = inputParser;
p.addRequired('A');
p.addOptional('y',[]);
p.addOptional('K',[]);
p.addOptional('N',[]);
p.addOptional('C',[]);
p.addOptional('tail', 'b', @(c)ismember( c, {'b','l','r'} ));
p.addOptional('perm',  0, @(x) isnumeric(x)&isscalar(x));  % total number of permutations
p.addOptional('repsz', 100, @(x) isnumeric(x)&isscalar(x)); % block size
p.addOptional('test', 'f', @(c)ismember( c, {'f','l'} )); %CR: Fisher or Liebermeister test
p.parse(varargin{:});
P = p.Results.C;
if isempty(p.Results.N)     % using fexact(X,y,options)
    X = p.Results.A;
    y = p.Results.y;
    if (~islogical(y) && ~all( ismember( y, [0 1])) ) 
        error('linstats:fexact:InvalidArgument', 'y must be in (0,1)' );
    end
    if (~islogical(X) && ~all( ismember( X(:), [0 1])) ) 
        error('linstats:fexact:InvalidArgument', 'X must be in (0,1)' );
    end

    y = logical(y);
    N = sum(y);           % number of cases
    M = length(y);
    K = sum(X,1)';
    x = X'*sparse(y);  % in unofficial testing using timeit. This was faster than sum(X(y==1,:))
    % and faster than sum(bsxfun(@eq,X,y))
else % using fexact( a,M,K,N ...)
    x = p.Results.A; 
    M = p.Results.y;
    N = p.Results.N;
    K = p.Results.K;
end

% force inputs to be column vectors so outputs are too
x = x(:);
K = K(:);

switch p.Results.tail
    case 'l'; tail = -1;
    case 'r'; tail = 1;
    case 'b'; tail = 2;
end

if p.Results.test == 'l' %CR added support for Liebermeister Test
    %Seneta E, Phipps MC (2001) On the Comparison of Two Observed Frequencies. Biometrical Journal, 43(1):23?43, 2001.
    %Phipps MC (2003)Inequalities between hypergeometric tails. Journal of Applied Mathematics and Decision Sciences. 7(3): 165-174
    %http://www.emis.de/journals/HOA/ADS/Volume7_3/174.pdf
    %p = fexact(5,31,15,6,'tail','r') %see Phipps p170-> 0.0721
    %p = fexact(5,31,15,6,'tail','r','test','l') %see Phipps p170-> 0.0345
    %p = fexact(10,31,15,25, 'tail','l')
    %p = fexact(10,31,15,25, 'tail','l','test','l')
    if tail == -1 %left tail: increment b,c
        %x = x+1;
        disp('left');
        N = N +1;
        K = K +1;
        M = M + 2;
    elseif tail == 1 %right tail: increment a,d
        x = x+1; 
        N = N +1;
        K = K +1;
        M = M + 2;
    else
        disp('Computing Fisher Exact values: Liebermeister measure only described for computing right tail.');
    end
end

if N==0 || N==M
    pval = ones(length(x),1);
    return;
end;

if isscalar(M) && isscalar(N) && isscalar(K) && isscalar(x) && ~(p.Results.perm > 1)
    pval = doTest( x, M, K, N, tail );
    pval(pval>1) = 1;   % fix roundoff errors
    return;
end

if isempty(P)
    uk = unique(K);
    P = getLookup( M,N,tail,uk);
end

pval = P( sub2ind(size(P), x+1, K+1));


if p.Results.perm > 1
    pperm   = doPerm(X,y,M,K,P,p.Results.perm,p.Results.repsz);
    [fx xx] = mecdf(pperm(:));
    xx(end) = max(pval);
    pval    = interp1( xx, fx, pval);
end

pval(pval>1) = 1;   % fix roundoff errors
%end fexact()

function pperm = doPerm(X,y,M,K,P,nperms,repsz)
% returns pperm a matrix with each element representing the smallest p
% value from a complete study where the rows of y have been randomly
% permuted
nperms = ceil(nperms/repsz); 
pperm = ones(repsz,nperms); % each element is from a permutation test. It is the minimum pvalue among all variables (columns of X)

% sub2ind was previously 85% of the execution time. 
% this method computes the index directly rather than
% computing x and then calliung sub2ind. It is much
% faster because it uses the fast multiple engine
% precompute K*nrows +1, then add x after it is calculated
X = [ K*size(P,1)+1 X' ]; 
B = ones(1,repsz);

for permi = 1:nperms
    % generate new random vector y, by shuffling the rows (do not change
    % the counts of M,K or M)
    [~, i] = sort( randn(M,repsz) ); % note randn was faster than rand when I checked
    pperm(:,permi) = min( P( X*[B;y(i)] ),[], 1); % find the smallest p-value among all separate tests for each permutation
end
%end pperm()

function P = doTest( x, M, N, K, tail)
minx = max(0,N+K-M);
maxx = min(N,K);
if x < minx
   error( 'x cannot be smaller than max(0,N+K-M)' );
end
if x > maxx;
   error( 'x cannot be larger than min(N,K)');
end
%log p(x|M,N,K) = Q(K,x)+Q(M-K,N-x)-Q(M,N) where Q(n,k) = lnchoosek(n,k)
if tail == -1
   a = minx:x;
   logP = lnchoosek( K, a) + lnchoosek(M-K,N-a) - lnchoosek(M,N);
   P = sum(exp(logP));
elseif tail == 1
   a = x:maxx;
   logP = lnchoosek( K, a) + lnchoosek(M-K,N-a) - lnchoosek(M,N);
   P = sum(exp(logP));
else
   a = minx:maxx;
   logP = lnchoosek( K, a) + lnchoosek(M-K,N-a) - lnchoosek(M,N);
   pcrit = logP(x-minx+1);
   P = exp(logP);
   k = P - exp(pcrit) <= M*length(P)*eps(max(P));
   P  = sum(P(k));
end

%end doTest()

function P = getLookup( M, N, tail, K)
% get Pvalue lookup table
%    returnns P(x+1,K+1) = log p(x,|M,N,K);
%

F = gammaln( 1:M+1);  % used to compute factorials. For small problems
% It is overkill to generate all factorials up
% to M, but it takes on .001 seconds to compute
% all of them to 1..10,001, so we can afford
% it. Memory-wise may be more of an issue, and
% a reason to change. 

% --- old way - compute all possible lookups for nchoosek
% Also we sometimes don't need a full lookup for ln
% nchoosek. we only lookup 
%     log p(x|M,N,K) = Q(K,x)+Q(M-K,N-x)-Q(M,N) 
%                where Q(n,k) = lnchoosek(n,k)
% n = repmat( n, size(k,1), 1);
% ind = bsxfun( @plus, n, (0:M-1)*M);
% [k, n, ind] = tri2sqind(M);
% L = zeros(M);
% L(ind) = lnchoosek( n(:), k(:), F);
% L = blkdiag(0,L');      % pad first row and col with 0s for n=0, k=0

% --- new way 
% build lookup table only for values we will need.
% We need Q(K,x)+Q(M-K,N-x)-Q(M,N) where Q(n,k) = lnchoosek(n,k)
% we also need all columns of Q for all K, M-K and M;
n = unique( [K(:); M-K(:); M] )'; % columns of initial (pre-transpsoed L);
k = repmat( (1:M)', 1, length(n)); % all rows
n = repmat( n, size(k,1), 1);

% L = zeros(M+1);           % pretty slow!
% L(ind) = lnchoosek( n(:), k(:), F);
% very much faster
L = sparse( n(:)+1, k(:)+1, lnchoosek( n(:), k(:), F) );


% to lookup log nchoosek(n,k) use L(n+1,k+1);
% create a lookup table for Log P
% to lookup log p(x,|M,N,K) use P(x+1,K+1);
% P initially contains 1s for valid values of x
% x is in [min( 0, N+K-M), ..., min(N,K)]
% min(N,K) is on the diagonal. the first
% column that contains a 0 in the first row will be when the minimum legal
% value of x is 1, so N+K-M == 1. K will equal M-N+1, which is the column
% M-N+2

%P will = log p(x|M,N,K) = Q(K,x)+Q(M-K,N-x)-Q(M,N) where Q(n,k) = lnchoosek(n,k)
% P = spdiags(ones(N+1,M+1),0:M-N,N+1,M+1);       % pretty slow! 
% but the whole point of this next section is just to get indices MKind,
% NXind, Kind ad Xind (P is overwritten below). There is likely a much
% better way
P = tril(triu(ones(N+1,M+1)),M-N);      % faster!
ki = ind2logical(n(1,:)'+1, size(P,2)); % only use the required lookup columns into L
P(:,~ki) = 0;                           
[Xind, Kind] = find(P);                 
MKind = M-Kind+2;
NXind = N-Xind+2;

q = sub2ind(size(P),Xind,Kind);
P = nan(size(P));

% FIXED: the next step is slow (function %time). I now use only the lookups
% required 
%   sub2ind 17%
%   L(k1)   45%
%   exp(sum(L)) 38%
% 

P(q) = exp( L(sub2ind(size(L),Kind,Xind))    + ... Q(K,x)
    L(sub2ind(size(L),MKind,NXind))  + ... Q(M-K, N-x)
    -L(M+1,N+1));                     % ... Q(M,N)
if tail == -1
    P = cumsum(tril(P, M-N));
elseif tail == 1
    P = flipud(triu(P));
    P = cumsum(P);
    R = N-(-1:size(P,1)-2)';
    P = P(R,:);
else
    % This function has been tested against two reference algorithm using all M
    % and N up to 30. The reference algorithms differ in how they find the
    % find the tail on the other side of the distribution.
    Pt = P(:,n(1,:)+1);
    [q, r] = size(Pt);
    [Pt, idx] = sort(Pt,1);     % sort order, called rank, will be transformed to rank below
    
    % FIXED: the next line is kind of slow (25% of all running time). Now I
    % only calculate results for the columns of P that will be used
    ties = [abs(diff(Pt)) < M*numel(P)*eps(Pt(1:end-1,:));
        zeros(1,r)];
    R = repmat( (1:q)',1,r);
    R(col2ind(idx)) = R+ties;
    Pt = cumsum(Pt,1);          % reuse P to store cumulative P
    P(:,n(1,:)+1) = Pt(col2ind(R));
end
%end getLookup()

function q = lnchoosek( n, k, f)
% lnchoosek natural log of nchoosek
%
% usage
%       q = lnchoosek( n, k );
% usage
%       q = lnchoosek( n, k, f); where f is the gammaln(1:max(n,k)+1);
q    = n-k;       % allocate space same size as n
i    = q>0;       % only compute for valid n and k

if nargin < 3
    % q is 0 when n = k because log nchoosek(a,a) = 1 for any a >= 0
    % otherwise it is given by ...
    q(i) = -log(q(i)) - betaln(k(i)+1,q(i));     % otherwise this is it
else
    q(i) =  f(n(i)+1) - f(k(i)+1) - f( q(i)+1 );
end

q(~i) = 0;
%end lnchoosek()

function i = col2ind(order,siz)
%col2ind converts column specific 1-based indices to element based indices
%i = col2ind(A, siz)
%   return i, a set of integer indices into a m x n matrix
%   A is a p x n set of column specific indices as you'd get from the second
%   output of sort
%   size is optional.
%        if present, P is set to SIZ(1), otherwise
%        P is set to the size of the first dimension of ORDER,
%
% example
%   [xsort order] = sort(x);
%   wrong = x(order);    % this isn't what you want
%   xsort = x(col2ind(order));   % this is

[m p] = size(order);
if nargin < 2
    q = m;
else
    q = siz(1);
end
oset = 0:q:(q*p-1);
i = order + repmat( oset, m, 1 );
%end col2ind()

function [i,j,k] = tri2sqind( m, k )
%TRI2SQIND subscript and linear indices for upper tri portion of matrix
%
% get indices into a square matrix for a vector representing a the upper
% triangular portion of a matrix such as those returned by pdist.
%
% [i,j,k] = tri2sqind( m, k )
%  If V is a hypothetical vector representing the upper triangular portion
%  of a matrix (not including the diagonal) and
%  M is the size of a square matrix and
%  K is an optional vector of indices into V then tri2sqind returns
%  (i,j) the subscripted indices into the equivalent square matrix.
%  K is an integer index into the equivalent square matrix
%
% Example
%  X = randn(5, 20);
%  Y = pdist(X, 'euclidean');
%  [i,j,k] = tri2sqind( 5 );
%  S = squareform(Y);
%  isequal( Y(:), S(k) );
%  Z = zeros(5);
%  Z(k) = Y;

max_k = m*(m-1)/2;

if ( nargin < 2 )
    k = (1:max_k)';
end;

if any( k > max_k )
    error('linstats:tri2sqind:InvalidArgument', 'ind2subl:Out of range subscript');
end;


i = floor(m+1/2-sqrt(m^2-m+1/4-2.*(k-1)));
j = k - (i-1).*(m-i/2)+i;
k = sub2ind( [m m], i, j );
%end tri2sqind()


function [fx x] = mecdf( y )
% MECDF my ECDF function to avoid dependendencies on stats toolbox
% assumes 0 <= y <= 1
% x always includes 0 and fx(x>1) = 1;
n = length(y);
[x i] = unique(sort(y));
fx = linspace( 0, 1, n+1 )';
fx = fx(i+1);
if x(1)>0
    x = [0;x];
    fx = [0;fx];
end

if x(end) < 1    % set up the maximum x value as Inf, so interp1 never fails 
    x  = [x;Inf];
    fx = [fx;1];
end
%end mecdf()



