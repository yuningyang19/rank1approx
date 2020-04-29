function [X,Y,Lam,k] = admm_rank1_d6(A,opt)
% solve min <A,X> s.t. ||X||=1, rank_{CP}(X)=1, X\in S^{n^4}
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
    factor = 1.5;
end
if isfield(opt,'X') == 1
    X = opt.X;
else
    X = zeros(size(A));
end
sz = size(X); n = sz(1);
Y = zeros(sz); Lam = -A;
opts.v0  = rand(opt.n^3,1);
count=0;

d = ndims(A);nd = 1:d;all_per = perms(nd);
flag = 0;
for k = 1: iter
    Lam = (Lam + permute(Lam,[4 5 6 1 2 3]))/2;
    Y0 = Y; Lam0 = Lam;
    Lam_tauY = -Lam0-tau*Y0; m_Lam_tauY = reshape(Lam_tauY,n^3,n^3); 
    m_Lam_tauY = (m_Lam_tauY+m_Lam_tauY')/2;
    [u,l] = eigs(m_Lam_tauY,1,'sa',opts);
    [lambda,min_i] = min(diag(l));
    v = u(:,min_i); opts.v0 = v;
    mX = v*v'; X = reshape(mX,sz); 
    t1 = A + Lam0 - tau*X;
    if abs(lambda)<= 1e-10
        flag = 1;
    end
    if flag == 0
        t2 = permute(t1 ,[1 2 4 3 5 6]); t3 = permute(t1 ,[1 2 4 5 3 6]); t4 = permute(t1 ,[1 2 4 5 6 3]);
        t5 = permute(t1 ,[1 4 2 3 5 6]); t6 = permute(t1 ,[1 4 2 5 3 6]); t7 = permute(t1 ,[1 4 2 5 6 3]);
        t8 = permute(t1 ,[4 1 2 3 5 6]); t9 = permute(t1 ,[4 1 2 5 3 6]); t10 = permute(t1 ,[4 1 2 5 6 3]);
        Y = -1/tau*(t1+t2+t3+t4+t5+t6+t7+t8+t9+t10)/10;         % Y = -1/tau(A+Lam-tauX);
    else
        Y = -1/tau* symmetrization(t1,all_per);
    end
    Lam = Lam0 - tau*factor*(X-Y);
    if   frob(X-Y)   <= eps
        break;
    end
    fprintf(1, repmat('\b',1,count));
    count = fprintf(1,'iter=%d, X-Y=%.4e',k,frob(X-Y));
end
fprintf('\n')
end
