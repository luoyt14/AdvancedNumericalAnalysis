function omega = bestOmega(A)
%寻找SOR迭代法最优的迭代因子w
%   input: H, output: omega

D = diag(diag(A));
L = D-tril(A);
U = D-triu(A);
step = 0.01;
rho = zeros(2/step-1,1);
k=1;
for w=step:step:2-step
    Lw=(D-w*L)\((1-w)*D+w*U);
    rho(k) = max(abs(eig(Lw)));
    k = k+1;
end

kbest = find(rho==min(rho));
omega = kbest(1)*step;

end

