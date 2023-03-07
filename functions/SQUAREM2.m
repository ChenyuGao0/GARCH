function paras = SQUAREM2(paras)
%%
n1 = paras.n1; n2 = paras.n2; Yn = paras.Yn; n = paras.n; n_max = max(n1,n2);

theta0 = [paras.omega0;paras.alpha0;paras.beta0;paras.sigma_t];
ct = paras.ct;
sigma_t = theta0(2+n1+n2:end);
ct0 = ct';
sigma0 = theta0(1) + ct0*theta0(2:1+n1+n2);
obj0 = sum(log(sigma_t))+...
    sum(Yn(n_max+1:end)./sigma_t)+paras.mu1*norm(sigma0-sigma_t)^2;

paras = obj4_s7(paras);
%debug_obj(paras,1);
paras = obj4_o1(paras);
%debug_obj(paras,1);
paras = obj4_ab1(paras);
%debug_obj(paras,1);
theta1 = [paras.omega0;paras.alpha0;paras.beta0;paras.sigma_t];
paras = obj4_s7(paras);
%debug_obj(paras,1);
paras = obj4_o1(paras);
%debug_obj(paras,1);
paras = obj4_ab1(paras);
%debug_obj(paras,1);
theta2 = [paras.omega0;paras.alpha0;paras.beta0;paras.sigma_t];
r = theta1 - theta0;
v = theta2 - theta1 - r;

alpha = min(norm(r,2)/(r'*v),-1);
alpha = max(alpha,-4);
while(true)
    theta3 = theta0 - 2*alpha*r + alpha^2*v;
    theta3 = projection(theta3,n1,n2);
    sigma_t = theta3(2+n1+n2:end);
    sigma_t1 = [paras.sigma_t0;sigma_t];

    for i = 1:n2
        ct(i+n1,:) = sigma_t1(n_max-i+1:n_max-i+n);
    end
    ct0 = ct';
    sigma0 = theta3(1) + ct0*theta3(2:1+n1+n2);
    obj1 = sum(log(sigma_t))+...
        sum(Yn(n_max+1:end)./sigma_t)+paras.mu1*norm(sigma0-sigma_t)^2;
    if(obj1<obj0)
        break;
    elseif(alpha==-1)
        break;
    else
        alpha = (alpha-1)/2;
    end
end
%fprintf('alpha = %3.8f\n', alpha);
paras.omega0 = theta3(1);
paras.alpha0 = theta3(2:1+n1);
paras.beta0 = theta3(2+n1:1+n1+n2);
paras.sigma_t = sigma_t;
paras.ct = ct;

%%
% paras = obj4_s8(paras);
% paras = obj4_o1(paras);
% paras = obj4_ab1(paras);
% theta4 = [paras.omega0;paras.alpha0;paras.beta0;paras.sigma_t];
% sigma_t = theta4(2+n1+n2:end);
% ct = paras.ct;
% ct0 = ct';
% sigma0 = theta4(1) + ct0*theta4(2:1+n1+n2);
% obj1 = sum(log(sigma_t))+...
%      sum(Yn(n_max+1:end)./sigma_t)+paras.mu1*norm(sigma0-sigma_t)^2;

%%
paras.diff = obj1 - obj0;
%fprintf('obj = %3.8f\n', obj1)
end

function theta_p = projection(theta,n1,n2)
eps = 1e-10;
omega0 = theta(1);
gamma0 = theta(2:1+n1+n2);
sigma_t = theta(2+n1+n2:end);
theta_p = zeros(size(theta));
theta_p(1) = max(eps,omega0);

l = -gamma0;
neg = l<0;
n_neg = sum(neg);
c = zeros(n1+n2,1);
if(sum(l(neg))+(1-eps)>=0)
    c(neg) = -l(neg);
else
    c_neg = c(neg);
    l_neg = l(neg);
    [lm,ind] = sort(l(neg),'descend');
    for i = 0:n_neg
        rho = -(sum(l_neg(ind(i+1:end)))+(1-eps))/(n_neg-i);
        if(prod(rho+lm(i+1:end)<0)==1)
            c_neg(ind(i+1:end)) = -(rho+l_neg(ind(i+1:end)));
            break;
        end
    end
    c(neg) = c_neg;
end
gamma0 = c;
theta_p(2:1+n1+n2) = gamma0;
theta_p(2+n1+n2:end) = max(sigma_t,eps);
end