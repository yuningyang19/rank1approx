function [X,Y,Lam,iter,eigvec,fval] = admm_rank1(A,opt)
% solve min <A,X> s.t. ||X||=1, rank_{CP}(X)=1, X\in S^{n^d}
% by nonconvex ADMM
% if max<A,X> needs to be computed, please use -A when calling this function instead.
% 
addpath ../PROPACK
addpath ../utils
    d = ndims(A);
    if d>=7 || d<=2
        error('order not support right now!\n')
    end
    switch d
        case 3
            [X,Y,Lam,iter] = admm_rank1_d3(A,opt); 
        case 4
            [X,Y,Lam,iter] = admm_rank1_d4(A,opt); 
        case 5
            [X,Y,Lam,iter] = admm_rank1_d5(A,opt); 
        case 6
            [X,Y,Lam,iter] = admm_rank1_d6(A,opt); 
    end
    sz = size(A); n = sz(1);
    mX = reshape(X,n,n^(d-1));
    [eigvec,~,~] = lansvd(mX,1,'L');
    kp_eigvec = eigvec;
    for i = 2 : d
        kp_eigvec = kron(eigvec,kp_eigvec);
    end
    fval = -A(:)'*kp_eigvec;     % because   A is in fact -A.
 
end