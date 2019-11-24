function [sX] = symmetrization(X,all_per)


    sX = 0;
    len_per = length(all_per);
    for i = 1:len_per
        sX = sX + permute(X,all_per(i,:));        
    end
    sX = sX/len_per;

end