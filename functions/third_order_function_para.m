function root_pr = third_order_function_para(para)
a = para(:, 1); b = para(:, 2); d = para(:, 4);
root_pr = zeros(length(a), 1);
t1 = - b.^3./(27*a.^3) - d./(2*a);
temp = d.^2./(4*a.^2) + (b.^3.*d)./(27*a.^4);
pos = temp>0;
if(prod(pos) == 1)
    delta1 = sqrt(temp);
    root_pr = -(b./(3*a)) + power(t1+delta1, 1/3) + power(abs(t1-delta1), 1/3);
else
    delta1 = sqrt(temp(pos));
    root_pr(pos) = -(b(pos)./(3*a(pos))) + power(t1(pos)+delta1, 1/3) + power(abs(t1(pos)-delta1), 1/3);

    neg = ~pos;
    delta2 = sqrt(abs(temp(neg)));
    root_pr(neg) = real(-(b(neg)./(3*a(neg))) + power(t1(neg)+delta2*1i, 1/3) + power(t1(neg)-delta2*1i, 1/3));
end
if(sum(real(root_pr)<=0)~=0)
    root_pr(real(root_pr)<=0) = 1e-6;
end

end