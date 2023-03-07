clear
addpath './functions'
addpath './data/real_data/'
returns = load('BAC_return.mat');
n1 = 1; n2 = 1; n_max = max(n1,n2);
k0 = 12;
scale = 50;
returns = struct2array(returns)*scale;
residuals = returns - mean(returns);
residuals2 = residuals.^2;
n = length(residuals);
emp_n = round(0.2*n);
n = n - emp_n;
Yn = residuals(1:n);
emp_Yn = residuals(n+1:end);
V0 = 0.1*var(Yn)*ones(n_max,1);
options = optimoptions(@fmincon,'Algorithm','interior-point');
n = n - n_max;
paras = struct;
EstMdl = estimate(garch(n2,n1),Yn,'Options',options);
o_m = EstMdl.Constant;
a_m = cell2mat(EstMdl.ARCH);
b_m = cell2mat(EstMdl.GARCH);
rng(k0+300)
paras.omega_ini = 0.1*var(Yn);
paras.alpha_ini = rand(n1,1);
paras.beta_ini = rand(n2,1);
if(sum(paras.alpha_ini)+sum(paras.beta_ini)>1)
    temp = sum(paras.alpha_ini)+sum(paras.beta_ini);
    paras.alpha_ini = 0.8*paras.alpha_ini/temp;
    paras.beta_ini = 0.8*paras.beta_ini/temp;
end

paras.omega0 = paras.omega_ini;
paras.alpha0 = paras.alpha_ini;
paras.beta0 = paras.beta_ini;

paras.sigma_t0 = V0;
sigma_t = [paras.sigma_t0;zeros(n,1)];
Yn2 = Yn.^2;
gamma0 = [paras.omega0;paras.alpha0;paras.beta0];
ct = [1; zeros(n1,1); zeros(n2,1)];

for i = 1:n
    ct(2:1+n1) = Yn2(i:i+n1-1);
    ct(2+n1:end) = sigma_t(i:i+n2-1);
    sigma_t(i+n_max) = ct'*gamma0;
end
ct = zeros(n1+n2,n);
for i = 1:n1
    ct(i,:) = Yn2(n_max-i+1:n_max-i+n);
end
for i = 1:n2
    ct(i+n1,:) = sigma_t(n_max-i+1:n_max-i+n);
end
paras.ct = ct;
paras.debug = 0;
paras.sigma_t = sigma_t(n_max+1:end);
paras.n = n;
paras.n1 = n1;
paras.n2 = n2;
paras.Yn = Yn2;
eps = 1e-6;
t = 1;
a_t = zeros(n1,1);
b_t = zeros(n2,1);
o_t = 0;
paras.obj = zeros(20000,1);
paras.t0 = t; paras.mu1 = 0;
%paras = debug_obj(paras,1,1);
%[mua,mub] = tryfirststep(paras);
sigma0 = zeros(n+n_max,1);
sigma0(1:n_max) = V0;
sigma_t = [paras.sigma_t0;paras.sigma_t];
total_t = zeros(20000,1);

ssize(1) = Inf;
while(t<Inf)
    tic
    a_t = paras.alpha0;
    b_t = paras.beta0;
    o_t = abs(paras.omega0);


    paras.mu1 = 3.9e10;


    paras = SQUAREM2(paras);
    %fprintf('t = %d\n', t);
    t = t + 1;
    total_t(t) = toc;


    ait(t,:) = a_t;
    bit(t,:) = b_t;
    oit(t,:) = o_t;
    if(abs(paras.diff)<eps)
        break;
    end

end
T = sum(total_t);
% plot(paras.sigma_t)
% hold on
% plot(Vn(n_max+1:end))
% hold off
performance_real_data(o_t/scale^2,a_t',b_t',residuals/scale,V0/scale,emp_n);
performance_real_data(o_m/scale^2,a_m,b_m,residuals/scale,V0/scale,emp_n);
% function [a,b] = tryfirststep(paras0)
% c = '*t^0.5';
% for i = 2:7
%     b = num2str(i-1);
%     for j = 1:3
%         a = num2str(1+3*(j-1));
%         paras = paras0;
%         mu0 = strcat(a,'e',b,c);
%         paras.mu0 = mu0;
%         for t = 1:100
%             paras.t0 = t;
%             paras.mu1 = eval(paras.mu0);
%             paras.mu2 = paras.mu1;
%             paras = obj4_s8(paras);
%             paras = obj4_o1(paras);
%             paras = obj4_ab1(paras);
%             paras = debug_obj(paras,1,1);
%         end
%         ssize = paras.obj(t) - paras.obj(t-1);
%         if(ssize<0)
%             a = 1+3*(j-1);
%             b = i-1;
%             return
%         end
%     end
% end
% end