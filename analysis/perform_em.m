function perform_em(intresults)


load(intresults)

pi_hat = pihat_k_num./repmat(sum(pihat_k_num,1),num_classes,1,1);
p_hati = zeros(1,num_classes);
for i = 1:num_classes
    p_hati(i) = nnz(Tinit(:,i))./num_tasks;
end

middle_prod = zeros(num_tasks,num_classes);
tic;
for m = 1:num_tasks
    if mod(m,10000)==0
        m
    end
    inner_prod = ones(num_players,num_classes);
    for i = 1:num_classes
        for j = 1:num_classes
            inner_prod(:,i) = inner_prod(:,i).*squeeze(pi_hat(i,j,:)).^(task_player_nmj{j}(m,:)');
        end
    end
    middle_prod(m,:) = prod(inner_prod(~isnan(sum(inner_prod,2)),:));% prod(inner_prod);
end
toc
outer_prod = prod(p_hati.*middle_prod);

min_votes = min(num_votes_pertask);
% num_votes_tot = sum(num_votes);
tot_votes_used = min_votes*num_tasks;
mu_0 = cat_votes_tot./sum(num_votes_pertask);
%task_votes_trim = cellfun(@(x) x(1:min_votes,:),task_votes,'UniformOutput',false);
currcol = 26;
task_votes_trim = cellfun(@(x) x(1:min_votes,currcol),task_votes,'UniformOutput',false);
pd_em(task_votes_trim,mu_0(currcol),w_0(currcol),prop_votes(:,currcol)>0.5)
