function paras = obj4_o1(paras)
n = paras.n;
eps = 1e-7;
gamma0 = [paras.alpha0; paras.beta0];
ct = paras.ct;
k0 = ct'*gamma0;
w0 = -1/n*sum(k0-paras.sigma_t);
if(w0>eps)
    paras.omega0 = w0;
else
    paras.omega0 = eps;
end
if(paras.debug)
    debug_obj(paras);
    fprintf('o = %1.8f\n', paras.omega0)
end
end