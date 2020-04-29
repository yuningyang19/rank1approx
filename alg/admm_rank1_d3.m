function [X,Y,Lam,k] = admm_rank1_d3(A,opt)
% solve min <A,X> s.t. ||X||=1, rank_{CP}(X)=1, X\in S^{n^3}
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
if isfield(opt,'n') == 1
    n = opt.n;
else
    error('parameter n must be specified!\n')
end
if isfield(opt,'X') == 1
    X = opt.X;
else
    X = zeros(size(A));
end
sz = size(X);
Y = zeros(sz); Lam = Y;
count=0;
opts.p0  = rand(opt.n,1);

for k = 1: iter    
    X0 = X;  Y0 = Y; Lam0 = Lam;    
    tmp1 = Lam0 + tau*Y0 + A ; tmp2 = permute(tmp1,[2 1 3]); tmp3 = permute(tmp1,[2 3 1]);
    X = 1/tau*(tmp1+tmp2+tmp3)/3;    
    Lam_tauX = Lam0 - tau*X;   m_Lam_tauX = reshape(Lam_tauX,n,n^2);  
    [u,l,v] = lansvd(m_Lam_tauX,1,'L',opts); opts.p0 = u;
    mY = -u*v'; Y = reshape(mY,n,n,n);
    Lam = Lam0 - tau*(X-Y);
     if frob(X-Y) <= eps
        break;
     end    
    fprintf(1, repmat('\b',1,count));
    count = fprintf(1,'iter=%d, X-Y=%.4e',k,frob(X-Y));        
end
fprintf('\n')
end
