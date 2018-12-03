function [x,k,rho] = Jacobi(A,b,x0,er1)
%Jacobi迭代法求解线性方程组Ax=b的解
%   input:A,b
%   output:x,k
x = x0;
er = 1;
k=0;

% for j=1:n-1                           				
%     [~,p]=max(abs(A(j:n,j)));
%     p=p+j-1;
%     if p>j
%         t1=A(j,:);A(j,:)=A(p,:);A(p,:)=t1;
%         t2=b(j,:);b(j,:)=b(p,:);b(p,:)=t2;
%     end 
% end

D = diag(diag(A));
L = tril(A)-D;
U = triu(A)-D;
J = -inv(D)*(L+U);
f = D\b;
rho = max(abs(eig(J)));
if rho<1
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
end

