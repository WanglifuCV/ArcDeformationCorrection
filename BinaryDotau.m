function [BinariedDotau]=BinaryDotau(Dotau)
[w,h]=size(Dotau);
BinariedDotau=zeros(w,h);
for i=1:w
    for j=1:h
        if (Dotau(i,j)>200)
      BinariedDotau(i,j)=255;
        else
            BinariedDotau(i,j)=0;
        end
    end
end

