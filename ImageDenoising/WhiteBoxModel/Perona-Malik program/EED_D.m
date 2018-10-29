function c=EED_D(gx, gy, type)
    # EED implementation
   grad2 = gx.*gx+gy.*gy;
   lambda1 = C(grad2);
   lambda2 = ones(size(grad2));
%   if type == 'x'
%   vector1 = gx ./ sqrt(grad2);
%   vector2 = gx ./ sqrt(grad2);
%   else 
%   vector1 = gy ./ sqrt(grad2);
%   vector2 = -gy ./ sqrt(grad2);
%   end
 
   vector1 = (gx + gy) ./ sqrt(grad2) ;
   vector2 = (gx - gy) ./ sqrt(grad2) ;
   
   vectorMat = [vector1, vector2];
   S = zeros(2*size(grad2));
   S(1:size(grad2)(1), 1:size(grad2)(1)) = lambda1;
   S(size(grad2)(1)+1:end, 1+size(grad2)(1):end) = lambda2;
   c = vectorMat * S * vectorMat';
end