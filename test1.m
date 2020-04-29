
% This test is to compute the largest eigenvalue of A.

clc; clear;
addpath ./funs
addpath ./PROPACK
d=4;  n=10;
sz = ones(1,d); sz = n*sz;

[A] = generate_Hilbert_tensor(n,d);



fa = frob(A);
A = A/fa;

opt.X = zeros(size(A)); opt.eps = 1e-4; opt.iter = 1000;   opt.n = n;
if mod(d,2) == 1
    opt.tau = 0.5;
else
    opt.tau = 0.1;
end

tic
[X,Y,Lam,iter,eigvec,fval] = admm_rank1(-A,opt);
if mod(d,2) == 1
    if fval < 0
        eigvec = -eigvec;
        fval = -fval;
    end
end
fval = fval*fa;     %rescalling
t1=toc;



if mod(d,2)==1
    mlam = reshape(Lam,n^(floor(d/2)) ,n^(ceil(d/2)));
    e = svd(mlam); maxeig = max(e)*fa;
else
    mlam = reshape(Lam,n^(d/2),n^(d/2));mlam=(mlam+mlam')/2;
    e = eig(mlam); maxeig = max(e)*fa;
end
if frob(fval-maxeig) <= 0.0001
    isglobal = 1;
else
    isglobal = 0;
end

fprintf('nonconvex: time = %.4f,  obj value = %.4f, iter = %d, isglobal = %d\n',t1,fval, iter,isglobal);







