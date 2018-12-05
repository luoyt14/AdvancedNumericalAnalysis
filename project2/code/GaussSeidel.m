function [x,k,rho,erlist] = GaussSeidel(A,b,x0,er1)
%GaussSeidel迭代法求解线性方程组Ax=b的解
%   input:A,b
%   output:x,k
x = x0;
er = 1;
erlist = [];
k=0;
xt = ones(size(x0));

D = diag(diag(A));
L = tril(A)-D;
U = triu(A)-D;
GS = -inv(D+L)*U;
f = (D+L)\b;
rho = max(abs(eig(GS)));
if rho<1+1e-6
    fprintf(1, 'GS迭代法收敛\n');
    while er>er1
        k = k+1;
        x1 = GS*x0+f;
        er = norm(x1-xt,2)/norm(xt,2);
        erlist = [erlist;er];
        x0 = x1;
    end
    x = x0;
else
    fprintf(1, 'GS迭代法不收敛\n');
end

% for m=1:20000
%     x1 = GS*x0+f;
%     x0 = x1;
% end
% x = x0;

end

