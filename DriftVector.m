function [F] = DriftVector(t,Y)

%   Purpose
%   =======
%   Return drift vector
%
%
%   IN
%   ==
%   1) t - time
%   2) Y - N-by-1 vector of solution at time t
%
%   OUT
%   ===
%   F - N-by-1 drift vector
%

%    F = [ Y(1); -Y(2) ];


    F = zeros(length(Y),1);

    % CLE
     c(1) = 1.0d3;
     c(2) = 1.0d3;
     c(3) = 1.0d-5;
     c(4) = 10.0d0;
     c(5) = 1.0d0;
     c(6) = 1.0d6;

%      stochiometric coefficients
     nu(1,1) = -1.d0;
     nu(1,2) = -1.d0;
     nu(1,3) =  1.d0;

     nu(2,1) =  1.d0;
     nu(2,2) =  1.d0;
     nu(2,3) = -1.d0;

     nu(3,1) = -1.d0;
     nu(3,2) =  1.d0;
     nu(3,3) = -1.d0;

     nu(4,1) =  1.d0;
     nu(4,2) = -1.d0;
     nu(4,3) =  1.d0;

     nu(5,1) =  1.d0;
     nu(5,2) = -1.d0;
     nu(5,3) = -1.d0;

     nu(6,1) = -1.d0;
     nu(6,2) =  1.d0;
     nu(6,3) =  1.d0;

     alpha(1) = c(1) * Y(1) * Y(2);
     alpha(2) = c(2) * Y(3);
     alpha(3) = c(3) * Y(1) * Y(3);
     alpha(4) = c(4) * Y(2);
     alpha(5) = c(5) * Y(2) * Y(3);
     alpha(6) = c(6) * Y(1);


   F(1) = nu(1,1)*alpha(1) + nu(2,1)*alpha(2) + nu(3,1)*alpha(3) + nu(4,1)*alpha(4) + nu(5,1)*alpha(5) + nu(6,1)*alpha(6);
   F(2) = nu(1,2)*alpha(1) + nu(2,2)*alpha(2) + nu(3,2)*alpha(3) + nu(4,2)*alpha(4) + nu(5,2)*alpha(5) + nu(6,2)*alpha(6);
   F(3) = nu(1,3)*alpha(1) + nu(2,3)*alpha(2) + nu(3,3)*alpha(3) + nu(4,3)*alpha(4) + nu(5,3)*alpha(5) + nu(6,3)*alpha(6);
   
end















