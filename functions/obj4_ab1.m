function paras = obj4_ab1(paras)
n1 = paras.n1;
n2 = paras.n2;
n_max = max(n1,n2);

sigma_t0 = paras.sigma_t0;
sigma_t = paras.sigma_t;
sigma_t = [sigma_t0;sigma_t];

eps = 1e-7;
gamma0 = [paras.alpha0; paras.beta0];

ct = paras.ct;

l1 = 2*ct*(paras.omega0-sigma_t(1+n_max:end));
k0 = norm(ct,'fro')^2;
L = ct*ct';

M = k0*eye(n1+n2);
k0 = 2*k0;
x1 = 2*(L-M)*gamma0;
l = x1+l1;

neg = l<0;
n_neg = sum(neg);
c = zeros(n1+n2,1);
if(sum(l(neg))+k0*(1-eps)>=0)
    c(neg) = -l(neg)/k0;
else
    c_neg = c(neg);
    l_neg = l(neg);
    [lm,ind] = sort(l(neg),'descend');
    for i = 0:n_neg
        rho = -(sum(l_neg(ind(i+1:end)))+k0*(1-eps))/(n_neg-i);
        if(prod(rho+lm(i+1:end)<0)==1)
            c_neg(ind(i+1:end)) = -(rho+l_neg(ind(i+1:end)))/k0;
            break;
        end
    end
    c(neg) = c_neg;
end
alpha0 = c(1:n1);
beta0 = c(n1+1:end);
paras.alpha0 = alpha0;
paras.beta0 = beta0;
% if(paras.debug)
%     sigma0 = zeros(n, 1);
%     for i = 1:n
%         st1 = i+(n_max-n1);
%         st2 = i+(n_max-n2);
%         ct = [flip(Yn(st1:st1+n1-1));flip(sigma_t(st2:st2+n2-1))];
%         gamma0 = [paras.alpha0;paras.beta0];
%         sigma0(i) = paras.omega0 + gamma0'*ct;
%     end
%     obj = sum(paras.sigma_t./paras.sigma_t1)-n+sum(log(paras.sigma_t1))+sum(Yn(2:end)./paras.sigma_t)+paras.mu1*norm(sigma0-paras.sigma_t)^2;
%     fprintf('obj = %f\n', obj)
% end
%fprintf('%1.8f, %1.8f\n', alpha0, beta0)
if(paras.debug)
    debug_obj(paras);
    fprintf('a = %1.8f, b = %1.8f\n', alpha0, beta0)
end
end