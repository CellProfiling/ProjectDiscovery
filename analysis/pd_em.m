function pd_em(task_votes,mu_0,w_0,xvec)

convergence = 0;
zeta = repmat(0.001,1,length(w_0));

%Y = expert currated HPA labels
num_tasks = size(task_votes,1);

%init weights
% mu_i0 = player_votes./numplayers;
% mu_i = p(yi=1|Yi,xi,theta)
% mu_i = ai*pi/(ai*pi+bi*(1-pi))
% p_y1 = sigmiod(wj'xi)
% mnrfit(

% coin_a%sensitivity - TPR - pr(yi=1|y=1)
% coin_b%specificity - (1-FPR) - pr(yi=0|y=0)
% gamma%threshold for decision - ROC= gamma from -inf to inf

%pr(y=1|x,w) = sigmoid(w'x) - logistic sigmoid
%set up iterative weights 
num_votes = size(task_votes{1},1);
num_classes = size(task_votes{1},2);
if size(mu_0,1)==1
    mu_n = repmat(mu_0,num_tasks,1);
else
    mu_n = mu_0;
end
% xvec = ones(1,num_tasks);
ai = zeros(num_tasks,num_classes);
bi = zeros(num_tasks,num_classes);
wt = repmat(w_0,1,num_votes);
% pi = sigmoid(wt'*xvec);


while convergence == 0
    
    
    
    %M-step
    alpha_num_j = zeros(num_votes,num_classes);
    beta_num_j = zeros(num_votes,num_classes);
    for i = 1:num_tasks
        num_votes = size(task_votes{i},1);
        alpha_num_j = alpha_num_j+repmat(mu_n(i,:),num_votes,1).*task_votes{i};
        beta_num_j = beta_num_j+repmat((1-mu_n(i,:)),num_votes,1).*(1-task_votes{i});
%         beta_j =
    end
    alpha_j = alpha_num_j./repmat(sum(mu_n,1),num_votes,1);
    beta_j = beta_num_j./repmat(sum((1-mu_n),1),num_votes,1);
    
    for i = 1:num_tasks
        ai(i,:) = prod((alpha_j.^task_votes{i}).*(1-alpha_j).^(1-task_votes{i}));
        bi(i,:) = prod((beta_j.^(1-task_votes{i})).*(1-beta_j).^(task_votes{i}));
    end
    
%     wt_1 = wt-zeta.*
    
%     curr_sigmoid = sigmoid(wt'*xvec);
    grad_w = zeros(num_votes,1);
    H_w = 0;
    %gradient step
    for i = 1:length(xvec)
        xi = xvec(i);
        grad_w = grad_w+(mu_n(i)-sigmoid(wt'*xi).*xi);
        H_w = H_w+sigmoid(wt'*xi).*(1-sigmoid(wt'*xi))*xi*xi';
    end
    
    %hessian
    H_w = sum(((curr_sigmoid).*(1-curr_sigmoid))*xvec*xvec');
    
    pi = sigmoid(w'*xvec);
    
    %E-step
    mu_n = (ai.*pi)/((ai*pi)+bi*(1-pi));
    
end

% go through each task and update weights 
% for task = 1:length(tasks)
%     annotators = tasks(task).annotators;
%     
%     
%     pi = sigmoid(w_i'xi);
% 
%    for ann = 1:length(annotators)
%        %take the product over the annotators
%        yij = annotators(ann);
%        alpha_j = sum(mu_i*yij)/sum(mu_i);
%        beta_j = sum((1-mu_i)*(1-yij))/sum(1-mu_i);
%        a_i = prod((alpha_j.^yij)*(1-alpha_j).^(1-yij))
%        b_i = prod((beta_j.^(1-yij))*(1-beta_j).^(yij))
%        
%    end
%    %take the product over the number of tasks
%     prD_theta = prod(a_i*pi+b_i*(1-pi));  
%    %theta ml = {alpha,beta,w} = argmax_theta(ln(prD_theta))
%    
%    e_lnpr = sum(mu_i*ln(pi)*a_i+(1-mu_i)*ln(1-pi)*b_i);
%    lnPr = sum(yi*ln(pi)*a_i+(1-yi)*ln(1-pi)*b_i);
%    mu_i = (a_i*pi)/(a_i*pi+b_i*(1-pi));
% 
%    
%    g = sum((mu_i-sigmoid(w'*xi))*xi);
%    H = -sum(sigmoid(w'*xi)*(1-sigmoid(w'xi))*xi*xi');
%    mu_new = mu_i-zeta*1/H*g; 
%    
% end



end
