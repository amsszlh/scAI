%% select high variable genes (HVGs) based on its average expression and Fano factor
function [id,mu,F_scaled] = HVGs(M,low_mu,high_mu,low_F,show)
if ~exist('low_mu','var') || isempty(low_mu)
    low_mu = 0.01; % select gene above this average expression level
end
if ~exist('high_mu','var') || isempty(high_mu)
    high_mu = 3.5; % select gene below this average expression level
end
if ~exist('low_F','var') || isempty(low_F)
    low_F = 0.5; % select gene above this Fano factor
end
if ~exist('show','var') || isempty(show)
    show = 0; 
end
% M is the raw count data matrix
% id is the set of select gene position

%% select HVGs: calculates the average expression and Fano factor for each gene, places these genes into bins,
% and then calculates a z-score for Fano factor within each bin.
mu = log(1+mean(exp(M)-1,2)); 
F = log(var(exp(M)-1,0,2)./mean(exp(M)-1,2));
mu(isnan(mu)) = 0;F(isnan(F)) = 0;
mu(isinf(mu)) = 0;F(isinf(F)) = 0;
[Y,E] = discretize(mu,20);
idx = setdiff(1:20,unique(Y));
if ~isempty(idx)
    E(idx+1) = [];
    Y = discretize(mu,E);
end
mean_y = grpstats(F,Y,'mean');
sd_y = grpstats(F,Y,'std');
F_scaled = (F - mean_y(Y))./sd_y(Y);
F_scaled(isnan(F_scaled)) = 0;
id = find((mu > low_mu) & (mu < high_mu) & F_scaled > low_F);
if show == 1
    figure
    scatter(mu,F_scaled,'k.')
    xlabel('Average expression');
    ylabel('Fano factor')
end