function [x,k,rho] = Jacobi(A,b,x0,er1)
%Jacobi迭代法求解线性方程组Ax=b的解
%   input:A,b
%   output:x,k
x = x0;
er = 1;
k=0;

D = diag(diag(A));
L = tril(A)-D;
U = triu(A)-D;
J = -inv(D)*(L+U);
f = D\b;
rho = max(abs(eig(J)));
if rho<1+1e-6
    fprintf(1, 'Jacobi迭代法收敛\n');
    while er>er1
        k = k+1;
        x1 = J*x0+f;
        er = norm(x1-x0,inf);
        x0 = x1;
    end
    x = x0;
else
    fprintf(1, 'Jacobi迭代法不收敛\n');
end

% for m=1:100
%     x1 = J*x0+f;
%     x0 = x1;
% end
% x = x0;
end

