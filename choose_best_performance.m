function best_one = choose_best_performance(result)
n = length(result); 
W_difference_corr = zeros(2,n);
for m = 1:n
    W1 = result{1,m}.W1; W2 = result{1,m}.W2;
    W_difference_corr(1,m) = X_difference_corr(W1);
    W_difference_corr(2,m) = X_difference_corr(W2);
end
W_difference_corr = sum(W_difference_corr)/2;
[~,Ws_index] = sort(W_difference_corr);
best_one = result{1,Ws_index(5)};
disp(['The best seed is ',num2str(Ws_index(5))])

% sum of correlation 
function X_difference_value = X_difference_corr(X)
%index = X_difference_matrix<0;X_difference_matrix(index) = X_difference_matrix(index)+1;
n = size(X,2);
if sum(X(:)) == 0
    X_difference_matrix = zeros(n,n);
else
    X_difference_matrix = corr(X);
end
X_difference_value = 0;
for i = 1:n
    for j = i+1:n
        X_difference_value = X_difference_value+X_difference_matrix(i,j);
    end
end