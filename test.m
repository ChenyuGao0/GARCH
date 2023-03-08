clear
addpath './functions'
load('./data/sim_data/data_norm400_(2,3).mat')
n1 = length(a); n2 = length(b);
k0 = 72;
Yn = Yn100(:,k0);
Vn = Vn100(:,k0);
n_max = max(n1,n2); n = length(Yn) - n_max;
paras = struct;
rng(k0+429+3)
paras.omega_ini = var(Yn)*0.1;
paras.alpha_ini = rand(n1,1);
paras.beta_ini = rand(n2,1);
if(sum(paras.alpha_ini)+sum(paras.beta_ini)>1)
    temp = sum(paras.alpha_ini)+sum(paras.beta_ini);
    paras.alpha_ini = 0.8*paras.alpha_ini/temp;
    paras.beta_ini = 0.8*paras.beta_ini/temp;
end
clearvars temp
paras.omega0 = paras.omega_ini;
paras.alpha0 = paras.alpha_ini;
paras.beta0 = paras.beta_ini;

paras.sigma_t0 = Vn(1:n_max);
sigma_t = [paras.sigma_t0;zeros(n,1)];
Yn = Yn.^2;
gamma0 = [paras.omega0;paras.alpha0;paras.beta0];
ct = [1; zeros(n1,1); zeros(n2,1)];

for i = 1:n
    ct(2:1+n1) = flip(Yn(i:i+n1-1));
    ct(2+n1:end) = flip(sigma_t(i:i+n2-1));
    sigma_t(i+n_max) = ct'*gamma0;
end
ct = zeros(n1+n2,n);
for i = 1:n1
    ct(i,:) = Yn(n_max-i+1:n_max-i+n);
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
paras.Yn = Yn;
eps = 10^(-6);
t = 1;
a_t = zeros(n1,1);
b_t = zeros(n2,1);
o_t = 0;
paras.obj = zeros(20000,1);
paras.t0 = t; paras.mu1 = 0;
%paras = debug_obj(paras,1,1);
%[mua,mub] = tryfirststep(paras);
sigma0 = zeros(n+n_max,1);
sigma0(1:n_max) = Vn(1:n_max);
sigma_t = [paras.sigma_t0;paras.sigma_t];
total_t = zeros(20000,1);
paras.q = 2;
paras.q0 = 0;
paras.R = [];
paras.V = [];
ssize(1) = Inf;
while(t<30000)
    tic
    a_t = paras.alpha0;
    b_t = paras.beta0;
    o_t = abs(paras.omega0);

     paras.mu1 = 1e5; %hyperparameter

    paras = SQUAREM2(paras);
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
plot(paras.sigma_t)
hold on
plot(Vn(n_max+1:end))
hold off
gamma0 = [o;a;b];