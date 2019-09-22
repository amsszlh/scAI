function [C,cluster,coph] = identify_consensus_cluster(result,K)
repeat = length(result);
Cs = cell(1,repeat);
for flag = 1:repeat
    h = result{1,flag}.H;
    H = h./repmat(sum(h),size(h,1),1);
    N = size(H,2); label = zeros(1,N);
    for j = 1:N
        sl = find(H(:,j) == max(H(:,j)));
        if ~isempty(sl)
            label(j) = sl(1);
        else
            label(j) = randi(K,1,1);
        end
    end
    % compute consensus matrix
    C = zeros(N);
    for j = 1:N
        for k = j:N
            if label(j) == label(k)
                C(j,k) = C(j,k)+1;
                C(k,j) = C(j,k);
            end
        end
    end
    Cs{1,flag} = C;
end
C = zeros(N);
for j = 1:repeat
    C = C+Cs{1,j};
end
C = C/repeat;
[~,cluster,~,coph] = nmforderconsensus0(C,K);