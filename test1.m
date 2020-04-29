
%   rng('default');
clc; clear;
addpath ./funs
addpath ./PROPACK
d=3;  n=20;
sz = ones(1,d); sz = n*sz;

[A] = generate_Hilbert_tensor(n,d);

% A =  randn(sz);
% d = ndims(A);nd = 1:d;all_per = perms(nd);
% A = symmetrization(A,all_per);

   
fa = frob(A); 
A = A/fa;

opt.X = zeros(size(A)); opt.eps = 1e-4; opt.iter = 1000;   opt.n = n;
if mod(d,2) == 1
    opt.tau = 0.5;
else
    opt.tau = 0.1;
end

tic
[X,Y,Lam,nonconvex_iter1,eigvec,fval] = admm_rank1(-A,opt,-fa);
if mod(d,2) == 1
    if fval < 0
        eigvec = -eigvec;
        fval = -fval;
    end
end
t1=toc;


abs(A(:)'*X(:)*fa-fval)     


% mlam = tens2mat(Lam,1:2,3:4);mlam=(mlam+mlam')/2;
% e = eig(mlam); maxeig = max(e);
% if frob(tmp-maxeig) <= 0.0001
%     isglobal = 1;
% else
%     isglobal = 0;
% end

fprintf('nonconvex time = %.4f,  value_org = %.4f, iter = %d\n',t1,fval, nonconvex_iter1);
fprintf('nie time = %.4f, value = %.4f\n',nie_time, nie_value);







