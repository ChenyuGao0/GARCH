function paras = obj4_s8_2(paras)
n = paras.n;
n1 = paras.n1;
n2 = paras.n2;
n_max = max(n1,n2);
mu1 = paras.mu1;
Yn = paras.Yn;

omega0 = paras.omega0;
alpha0 = paras.alpha0;
beta0 = paras.beta0;
sigma_t = paras.sigma_t;
sigma_t0 = paras.sigma_t0;
sigma_t1 = [sigma_t0; sigma_t];
x0 = -Yn(n_max+1:end);
dt = [-1; beta0];
l_max = dt'*dt;
M = l_max*eye(n2+1);
l_max = mu1*l_max;
D = dt*dt';
DM = D - M;
ct = zeros(n1+n2,n);
bt = [sigma_t1(n_max+1:end)';zeros(n2,n)];
b3 = [ones(1,n);zeros(n2,n)];
for i = 1:n1
    ct(i,:) = Yn(n_max-i+1:n_max-i+n);
end
for i = 1:n2
    ct(i+n1,:) = sigma_t1(n_max-i+1:n_max-i+n);
    bt(1+i,:) = ct(i+n1,:);
end
c0 = omega0 + ct(1:n1,:)'*alpha0;
b2 = DM*bt;
c1 = [c0';zeros(n2,n)];
for i = 1:n2
    b2(i+1,:) = [b2(i+1,i+1:end),zeros(1,i)];
    b3(i+1,1:n-i) = ones(1,n-i);
    c1(i+1,1:n-i) = c0(i+1:end)';
end
b2_ = sum(b2);

% sigma = c0 + ct(n1+1:end,:)'*beta0;
x2 = mu1*2*(b2_'+(dt'*c1)');
x3 = 2*l_max*sum(b3);

root_pr = zeros(n,1);
for i = 1:100
    root = roots([x3(i),x2(i),1,x0(i)]);
    root = root(root == real(root));
    p_root = root(root>0);
    if(length(p_root)>1)
        pos1 = p_root(1); pos2 = p_root(2);
        value_pr1 = log(pos1)-x0(i)./pos1+x2(i).*pos1+x3(i).*pos1.^2;
        value_pr2 = log(pos2)-x0(i)./pos2+x2(i).*pos2+x3(i).*pos2.^2;
        if(value_pr1>value_pr2)
            root_pr(i) = pos1;
        else
            root_pr(i) = pos2;
        end
    else
        root_pr(i) = p_root;
    end
end

paras.sigma_t = root_pr;
end