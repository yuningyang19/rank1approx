function [A] = generate_Hilbert_tensor(n,d)

   t=[];
   for i = 1: d
       t = [t n];
   end    
   A = zeros(t);
    
   if d == 3
       for i1= 1: n
           for i2 = 1: n
               for i3 = 1: n
                   A(i1,i2,i3) = 1/(i1+i2+i3-d+1);
               end
           end
       end
   elseif d == 4
       for i1= 1: n
           for i2 = 1: n
               for i3 = 1: n
                   for i4 = 1: n
                    A(i1,i2,i3,i4) = 1/(i1+i2+i3+i4-d+1);
                   end
               end
           end
       end
   elseif d == 5
      for i1= 1: n
           for i2 = 1: n
               for i3 = 1: n
                   for i4 = 1: n
                       for i5 = 1: n
                            A(i1,i2,i3,i4,i5) = 1/(i1+i2+i3+i4+i5-d+1);
                       end
                   end
               end
           end
       end       
   elseif d == 6
            for i1= 1: n
           for i2 = 1: n
               for i3 = 1: n
                   for i4 = 1: n
                       for i5 = 1: n
                           for i6 = 1: n
                                A(i1,i2,i3,i4,i5,i6) = 1/(i1+i2+i3+i4+i5+i6-d+1);
                           end
                       end
                   end
               end
           end
       end 
   end
%     for i = 1: r
%         b(:,i) = b(:,i)/norm(b(:,i));
%   
end