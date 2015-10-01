function [Y,Wiener] = Milstein(DriftVector,DiffusionMatrix,T,Y0,Wiener)

%   Purpose
%   =======
%   Find solution of the system of Ito stochastic equations with 
%   multi-channel non-commutative noise:
%
%      / Y1 \   / f1 \        / g11  g12 ... g1m \   / dW1 \
%      | Y2 |   | f2 |        | g12  g22 ... g2m |   | dW2 |
%      | .  | = | .  | * dt + |  .      .     .  | * |  .  |
%      | .  |   | .  |        |  .         .  .  |   |  .  |
%      \ Yn /   \ fn /        \ gn1  gn2 ... gnm /   \ dWm /
%
%       Yi(t0) = Yi0, i = 1..n
%
%
%   Method
%   ======
%   Milstein method on uniform time grid:
%                                                   / I1 \   
%                                                   | I2 |
%   yi[k+1] = yi[k] + h*fi[k] + [ gi1 gi2 ... gim ] | .  | + Tr( Jgi * A )   
%                                                   | .  |
%                                                   \ Im /
%   where
%
%        / dgi1/dx1  dgi1/dx2  ...  dgi1/dxn \
%        | dgi2/dx1  dgi2/dx2  ...  dgi2/dxn |
%  Jgi = |   .                  .       .    |
%        |   .                    .     .    |
%        \ dgim/dx1  dgim/dx2  ...  dgim/dxn /
%
%        / g11  g12 ... g1m \  / I11  I12 ... I1m \ 
%        | g21  g22 ... g2m |  | I21  I22 ... I2m | 
%   A =  |  .      .     .  |  |  .      .     .  | 
%        |  .         .  .  |  |  .         .  .  | 
%        \ gn1  gn2 ... gnm /  \ Im1  Im2 ... Imm / 
%
%   and Tr is a trace operator
%
%
%   IN
%   ==
%   1) DriftVector     - function handle that evaluates drift vector
%   2) DiffusionMatrix - function handle that evaluates matrix of 
%                        diffusion coefficients
%   3) tspan - vector of time points
%   4) Y0 - n-dimensional column vector with initial data
%   5) varargin - optional array of driving Wiener processes 
%                 (same as in the output below)
%
%
%   OUT
%   ===
%   T - K-dimensional column vector of time points 
%   Y - K-by-n solution array. Each row in Y corresponds to the solution 
%       at a time returned in the corresponding row of T
%   Wiener - M-by-K-dimensional array of the driving Wiener processes. 


    % number of equations
    N = max(length(Y0));
    
    % dimension of the noise
    M = size(DiffusionMatrix(1,Y0),2);
    
    % number of points in time discretization
    K = max(length(T));

    % step size
    dt = T(2) - T(1);

    % initialize solution array
    Y = zeros(N,K);

    Y(:,1) = Y0(:);
    
    
    % generate array of driving Wiener processes
    if ( nargin == 4 )
        Wiener = BrownianMotion(dt,K,M);
    end
    
    % loop in time
    for i = 2:K
        % generate vector of noise increments
        dW = Wiener(:,i) - Wiener(:,i-1);

        F = DriftVector(T(i-1),Y(:,i-1));
        G = DiffusionMatrix(T(i-1),Y(:,i-1));
        
        Ito = MultIto(dt,M,dW);
                
        % update solution
        Y(:,i) = Y(:,i-1) + F*dt + G*dW + MultItoPart2(T(i-1),Y(:,i-1));
    end
    
    
    function result = MultItoPart(t,X)
        result = zeros(N,1);
        GPrime = DiffusionJacob(t,X);
        B = G * Ito;
        for jj = 1:N
            bufGG = B .* GPrime(:,:,jj);
            result(jj) = sum(bufGG(:));
        end
    end

    function result = MultItoPart2(t,X)
        dx = 1e-6;
        result = zeros(N,1);
        B = (G * Ito)';
        newX = X;
        for jj = 1:N
            newX(jj) = X(jj) + dx;
            GPrime = ( DiffusionMatrix(t,newX) - G ) ./ dx;
            result = result + GPrime * B(:,jj);
            newX(jj) = X(jj);
        end
    end

    function result = MultItoPart3(t,X,ItoMatrix)
        dx = 1e-6;
        result = zeros(N,1);
        bufG = DiffusionMatrix(t,X);
        newX = X;
        for jj1 = 1:M
            for jj2 = 1:M
                for kk = 1:N
                    newX(kk) = newX(kk) + dx;
                    bufGPrime = (DiffusionMatrix(t,newX) - DiffusionMatrix(t,X))./dx;
                    result = result + ItoMatrix(jj1,jj2)*bufG(kk,jj1)*bufGPrime(:,kk);
                    newX(kk) = newX(kk) - dx;
                end
            end
        end
    end
    
    function J = Jacob()
        dx = 1e-6;
        J = zeros(N,N);
        X1 = Y(:,i);
        for ii = 1:N
            X1(ii) = Y(ii,i) + dx;
            J(:,ii) =  DriftVector(T(i),X1) - oldDrift;
            X1(ii) = Y(ii,i);
        end
        J = J / dx;
    end
    
    function Jg = DiffusionJacob(t,X)
    %   Purpose
    %   =======
    %   Find jacobians of rows of diffusion matrix:
    %   
    %               / dGi1/dX1  dGi2/dX1  ...  dGim/dX1 \
    %               | dGi1/dX2  dGi2/dX2  ...  dGim/dX2 |
    %   Jg(i,:,:) = |   .                   .       .   |
    %               |   .                     .     .   |
    %               \ dGi1/dXn  dGi2/dXn  ...  dGim/dXn /
    %
    %
    %   IN
    %   ==
    %   1) t - time
    %   2) X - N-by-1 vector of solution at time t
    %
    %   OUT
    %   ===
    %   Jg - N-by-M-by-N array. Each n-th N-by-M slice of the array
    %        corresponds to the transposed Jacobian of the n-th row 
    %        of the diffusion matrix
    %   
        dx = 1e-6;
        Jg = zeros(N,M,N);
        X1 = X;
        for jj = 1:N
            X1(jj) = X(jj) + dx;
            buf = DiffusionMatrix(t,X1) - DiffusionMatrix(t,X); 
            for ii = 1:N
                Jg(jj,:,ii) = buf(ii,:);
            end
            X1(jj) = X(jj);
        end
        Jg = Jg / dx; 
    end
end


























