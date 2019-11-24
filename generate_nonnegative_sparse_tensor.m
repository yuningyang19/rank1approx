function [A] = generate_nonnegative_sparse_tensor(n,d,sparse_level)

   t=[];
   for i = 1: d
       t = [t n];
   end
   A = tensor(rand(t));  A = double(symmetrize(A)); 
   
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