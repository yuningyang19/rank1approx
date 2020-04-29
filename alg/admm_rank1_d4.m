function [X,Y,Lam,k] = admm_rank1_d4(A,opt)
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
Y = zeros(sz);
Lam = -A;
opts.v0  = rand(opt.n^2,1); % opts.v0 = opts.v0/norm(opts.v0);
count=0;

d = ndims(A);nd = 1:d;all_per = perms(nd);
flag = 0;

for k = 1: iter    
    Y0 = Y; Lam0 = Lam;    
    Lam_tauY = -Lam0-tau*Y0;   m_Lam_tauY = reshape(Lam_tauY,n^2,n^2);    
    m_Lam_tauY = (m_Lam_tauY+m_Lam_tauY')/2;
    [u,lambda] = eigs(m_Lam_tauY,1,'sa',opts);
    v=u; opts.v0 = v;
    mX = v*v';    
    X = reshape(mX,sz);
    t1 = A + Lam0 - tau*X;
    if abs(lambda)<= 1e-10  
        flag = 1;
    end
    if flag == 0
        t2 = permute(t1 ,[1 3 2 4]); t3 = permute(t1 ,[1 3 4 2]);
        Y  = -1/tau*(t1  +t2+t3   )/3;
    else
        Y = -1/tau* symmetrization(t1,all_per);
    end
    Lam = Lam0 - tau*factor*(X-Y);    
    fprintf(1, repmat('\b',1,count));
    count = fprintf(1,'iter=%d, X-Y=%.4e',k,frob(X-Y));
    if   frob(X-Y)   <= eps
        break;
    end    
end
fprintf('\n')
end
