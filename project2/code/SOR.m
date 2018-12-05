function [x,k,rho,erlist] = SOR(A,b,x0,er1,w)
%SOR迭代法求解线性方程组Ax=b的解
%   input:A,b,x0,er1,w
%   output:x,k,rho
x = x0;
er = 1;
erlist = [];
k=0;
xt = ones(size(x0));

D = diag(diag(A));
L = D-tril(A);
U = D-triu(A);
Lw = (D-w*L)\((1-w)*D+w*U);
f = (D-w*L)\b*w;

rho = max(abs(eig(Lw)));
if rho<1+1e-6
    fprintf(1, 'SOR迭代法收敛\n');
    while er>er1
        k = k+1;
        x1 = Lw*x0+f;
        er = norm(x1-xt,2)/norm(xt,2);
        erlist = [erlist;er];
        x0 = x1;
    end
    x = x0;
else
    fprintf(1, 'SOR迭代法不收敛\n');
end

% for m=1:12441
%     x1 = Lw*x0+f;
%     x0 = x1;
% end
% x = x0;
end

