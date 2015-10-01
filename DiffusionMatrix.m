function G = DiffusionMatrix(t,Y)

%   Purpose
%   =======
%   Return diffusion matrix
%
%
%   IN
%   ==
%   1) t - time
%   2) Y - N-dimensional vector of solution at time t
%
%   OUT
%   ===
%   G - N-by-M matrix
%

%     G = [ Y(1),      0; ...
%              0,  2*Y(2) ];


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

     G(1,1) = nu(1,1) * sqrt(alpha(1));
     G(1,2) = nu(2,1) * sqrt(alpha(2));
     G(1,3) = nu(3,1) * sqrt(alpha(3));
     G(1,4) = nu(4,1) * sqrt(alpha(4));
     G(1,5) = nu(5,1) * sqrt(alpha(5));
     G(1,6) = nu(6,1) * sqrt(alpha(6));

     G(2,1) = nu(1,2) * sqrt(alpha(1));
     G(2,2) = nu(2,2) * sqrt(alpha(2));
     G(2,3) = nu(3,2) * sqrt(alpha(3));
     G(2,4) = nu(4,2) * sqrt(alpha(4));
     G(2,5) = nu(5,2) * sqrt(alpha(5));
     G(2,6) = nu(6,2) * sqrt(alpha(6));

     G(3,1) = nu(1,3) * sqrt(alpha(1));
     G(3,2) = nu(2,3) * sqrt(alpha(2));
     G(3,3) = nu(3,3) * sqrt(alpha(3));
     G(3,4) = nu(4,3) * sqrt(alpha(4));
     G(3,5) = nu(5,3) * sqrt(alpha(5));
     G(3,6) = nu(6,3) * sqrt(alpha(6));

end



















