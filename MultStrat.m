function StratInt = MultStrat(dt,M,ksi)

%   Purpose
%   =======
%   Find values of the multiple Stratonovich integrals of the form:
%                   /t  /s    j1   j2
%        J(j1,j2) = |   |   dW   dW
%                   /0  /0    
%
%   Method
%   ======
%   Karhunen-Loeve (Fourier) series expansion:
%       Ref - P.Kloeden "Numerical solution of stochastic differential
%                        equation", Chapter 5.8, Chapter 10.3
%
%   Additionally we utilize the property:
%
%       J(j1,j2) + J(j2,j1) = J(j1)*J(j2)
%           Ref - P.Kloeden "Numerical solution of stochastic differential
%                            equation", (10.3.15)
%
%
%   IN
%   ==
%   1) dt  - integrating time step
%   2) M   - dimension of the white noise
%
%   OUT
%   ===
%   Ito - M-by-M matrix with multiple Ito integrals
%
    
    sqrt2 = sqrt(2);

    % initialize matrix
    StratInt = zeros(M,M);

    % find number of partial sums required to get 1-st order convergence
    % Ref - Kloeden, (10.3.10)
    p = ceil(1/dt);  
    
    % additional coefficient
    ro = 0;
    for r = p:-1:1
        ro = ro + 1/(r*r);
    end
    ro = sqrt(1/12 - 0.5*ro/(pi*pi));
        
    ksi = ksi ./ sqrt(dt);
    mu  = randn(M,1);
    
    % non-diagonal entries
    for r = p:-1:1
        nu   = randn(M,1);
        zeta = randn(M,1);
        for i = 1:M
            StratInt(:,i) = StratInt(:,i) + ...
                            1/r * ( zeta(i)*( sqrt2*ksi    + nu    ) - ...
                                    zeta   *( sqrt2*ksi(i) + nu(i) ) );
        end
    end
    for i = 1:M
        StratInt(:,i) = StratInt(:,i) / pi + ... 
                        ksi(i)*ksi + 2*ro*( mu(i)*ksi - mu*ksi(i) );
    end
    StratInt = dt * 0.5 * StratInt;

end



























