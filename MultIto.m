function ItoInt = MultIto(dt,M,ksi)

%   Purpose
%   =======
%   Find values of the multiple Ito integrals of the form:
%                   /t  /s    j1   j2
%        I(j1,j2) = |   |   dW   dW
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
%       I(j1,j2) + I(j2,j1) = I(j1)*I(j2)
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
    
    ItoInt = MultStrat(dt,M,ksi);
    for i = 1:M
        ItoInt(i,i) = ItoInt(i,i) - 0.5 * dt;
    end
    
end



























