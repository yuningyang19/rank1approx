function [X,Y,Lam,k] = admm_rank1_d5(A,opt)
% solve min <A,X> s.t. ||X||=1, rank_{CP}(X)=1, X\in S^{n^5}
% L(X,Y,Lam) = <A,X> - <Lam, X-Y> + tau/2*||X-Y||^2
% Ideally, the limit satisfies X*=Y^*, S(Lam*)=A, X^*=argmin<Lam^*-tau*X^*,X>
if isfield(opt,'iter') == 1
    iter = opt.iter;
else
    iter = 1000;
end
if isfield(opt,'eps') == 1
    eps = opt.eps;
else
    eps = 1e-4;
end
tau = opt.tau;
if isfield(opt,'factor') == 1
    factor=opt.factor;
else
    factor = 1.;
end
if isfield(opt,'X') == 1
    X = opt.X;
else
    X = zeros(size(A));
end
sz = size(X); n = sz(1);
Y = zeros(sz); Lam = Y;
count=0;
opts.p0  = rand(opt.n^2,1);

d = ndims(A);nd = 1:d;all_per = perms(nd); flag = 0;
for k = 1: iter    
    X0 = X;  Y0 = Y; Lam0 = Lam;    
    tmp1 = Lam0 + tau*Y0 + A;    
    if flag == 0
        tmp2 = permute(tmp1,[1 3 2 4 5]); tmp3 = permute(tmp1,[1 3 4 2 5]); tmp4 = permute(tmp1,[1 3 4 5 2]);
        tmp5 = permute(tmp1,[3 1 2 4 5]); tmp6 = permute(tmp1,[3 1 4 2 5]); tmp7 = permute(tmp1,[3 1 4 5 2]);
        tmp8 = permute(tmp1,[3 4 1 2 5]); tmp9 = permute(tmp1,[3 4 1 5 2]); tmp10 = permute(tmp1,[3 4 5 1 2]);
        X = 1/tau* (tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8+tmp9+tmp10)/10 ;
    else
        X = 1/tau* symmetrization(tmp1,all_per);
    end    
    Lam_tauX = Lam0 - tau*X;   m_Lam_tauX = reshape(Lam_tauX,n^2,n^3); 
    [u,lambda,v] = lansvd(m_Lam_tauX,1,'L',opts); opts.p0 = u;
    mY = -u*v'; Y = reshape(mY,sz); 
    if abs(lambda)<= 1e-10 
        flag = 1;
    end
    Lam = Lam0 - tau*(X-Y);    
    if    frob(X-Y)    <= eps
        break;
    end    
    fprintf(1, repmat('\b',1,count));
    count = fprintf(1,'iter=%d, X-Y=%.4e',k,frob(X-Y));
end
fprintf('\n')
end
