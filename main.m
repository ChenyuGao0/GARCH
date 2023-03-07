clear
addpath './functions'
load(['./data/sim_data/' ...
    'data_norm300_(2,3).mat'])
n1 = length(a); n2 = length(b); n_max = max(n1, n2);
[n,N] = size(Yn100);
n = n - n_max;
as2 = zeros(N, n1);
bs2 = zeros(N, n2);
os2 = zeros(N);
res = zeros(N, n1+n2+1);
ts = zeros(N);
T = zeros(N);

for k0 = 1:N
    Yn = Yn100(:,k0);
    Vn = Vn100(:,k0);
    paras = struct;
    rng(k0+429+dice)
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

    paras.sigma_t0 = Vn(1:n_max);
    sigma_t = [paras.sigma_t0;zeros(n,1)];
    Yn = Yn.^2;
    gamma0 = [paras.omega0;paras.alpha0;paras.beta0];
    ct = [1; zeros(n1,1); zeros(n2,1)];

    for i = 1:n
        ct(2:1+n1) = Yn(i:i+n1-1);
        ct(2+n1:end) = sigma_t(i:i+n2-1);
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
    sigma_t_db = zeros(n, 100);
    sigma_t_db(:,1) = paras.sigma_t;
    eps = 1e-5;
    t = 1;
    total_t = zeros(20000,1);
    a_t = zeros(n1,1);
    b_t = zeros(n2,1);
    o_t = 0;
    paras.obj = zeros(20000,1);
    paras.t0 = t; paras.mu1 = 0;
    sigma0 = zeros(n+n_max,1);
    sigma0(1:n_max) = Vn(1:n_max);
    sigma_t = [paras.sigma_t0;paras.sigma_t];
    while(t<10000)
        tic
        a_t = paras.alpha0;
        b_t = paras.beta0;
        o_t = abs(paras.omega0);

        paras.mu1 = 1e4*1.0001^t;
        paras = SQUAREM2(paras);
        t = t + 1;
        total_t(t) = toc;

        if(abs(paras.diff)<eps)
            break;
        end

    end
    ts(k0) = t;
    os2(k0) = o_t;
    as2(k0,:) = a_t;
    bs2(k0,:) = b_t;
    res(k0,:) = [o_t,a_t',b_t'];
    T(k0) = sum(total_t);
    fprintf('k0 = %d\n', k0)
end
form = '.mat';
TimeNow = datestr(now,'mm-dd-HH-MM-SS');
obj0 = './result/SQUAREM_norm300_(2,3)';
filename = strcat(obj0,'_',TimeNow,form);
save(filename,'res','ts','total_t',"T");
