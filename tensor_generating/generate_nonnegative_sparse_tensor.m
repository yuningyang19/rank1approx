function [A] = generate_nonnegative_sparse_tensor(n,d,sparse_level)

t = ones(1,d); t = t*n;
A =   rand(t);
all_per = perms(1:d); A = symmetrization(A,all_per);

is_true = 1;
mu = 0.5;
while is_true
b=sum(A(:)>mu);
if b/n^d > sparse_level
mu = mu + 0.001;
else
is_true = 0;
end
end
A = A.* double(A>mu);



end
