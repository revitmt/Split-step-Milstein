function [T,Y,err] = Milstein(DriftVector,DiffusionMatrix,tspan,K,Y0)

%   Purpose
%   =======
%   Find solution of the system of Ito stochastic equations with 
%   m-dimensional multi-channel noise:
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
%   Milstein method:
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
%   1) DriftVector     - pointer to n-dimensional column vector of drift 
%                        coefficients
%   2) DiffusionMatrix - pointer to n-by-m matrix of diffusion coefficients
%   3) tspan = [T0; TN] - vector with interval of integration
%   4) K - number of intervals of integration
%   5) Y0 - n-dimensional column vector with initial data
%
%
%   OUT
%   ===
%   T - K-dimensional column vector of time points 
%   Y - K-by-n solution array. Each row in Y corresponds to the solution 
%       at a time returned in the corresponding row of T


    % number of equations
    N = max(length(Y0));
    
    % dimension of the noise
    M = size(DiffusionMatrix(1,Y0),2);

    % step size
    dt = (tspan(2)-tspan(1))/K;

    % initialize vector of time points
    T = zeros(1,K+1);

    % initialize solution array
    Y = zeros(N,K+1);
    
    T(1) = tspan(1);
    Y(:,1) = Y0(:);
    
    sqrt_dt = sqrt(dt);
    
    WinenerInt = zeros(M,1);

    time1 = 0;
    time2 = 0;
    
    % loop in time
    for i = 2:(K+1)
        % update time
        T(i) = T(1) + (i-1)*dt;
        
        % generate vector of the noise increment
        Wiener = randn(M,1);
        
        WinenerInt = WinenerInt + sqrt_dt*Wiener;

        G = DiffusionMatrix(T(i),Y(:,i-1));
        Ito = MultStrat(dt,M,Wiener);
        Ito(logical(eye(M))) = Ito(logical(eye(M))) - 0.5 * dt;
                
        % update solution
        Y(:,i) = Y(:,i-1) + dt*DriftVector(T(i),Y(:,i-1)) + sqrt_dt*G*Wiener + MultItoPart(T(i),Y(:,i-1),Ito);
        
       
%         if (min(Y(:,i)) < 0)
%             disp('negative value') ;
            for j=1:N
                if(Y(j,i)<0)
%                     Y(j,i) = 0;
%                     disp('negative value') ;
                end
            end
%         end
    end

    
    analitSol1 = exp(-2*T(K+1) + WinenerInt(1) - WinenerInt(2)) * cos(WinenerInt(3));
    analitSol2 = exp(-2*T(K+1) + WinenerInt(1) - WinenerInt(2)) * sin(WinenerInt(3));
    err(1) = (Y(1,K+1) - analitSol1);
    err(2) = (Y(2,K+1) - analitSol2);
    
    function result = MultItoPart(t,X,ItoMatrix)
        result = zeros(N,1);
        bufGPrime = DiffusionJacob(t,X);
        bufG = DiffusionMatrix(t,X)*ItoMatrix;
        for jj = 1:N
            bufGG = bufG.*bufGPrime(:,:,jj);
            result(jj) = sum(bufGG(:));
        end
    end

    function result = MultItoPart2(t,X,ItoMatrix)
        dx = 1e-6;
        result = zeros(N,1);
        bufG = DiffusionMatrix(t,X)*ItoMatrix;
        newX = X;
        for jj = 1:N
            newX(jj) = newX(jj) + dx;
            bufGPrime = (DiffusionMatrix(t,newX) - DiffusionMatrix(t,X))./dx;
            result = result + bufGPrime*bufG(jj,:)';
            newX(jj) = newX(jj) - dx;
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
    %   Find jacobians of the rows of diffusion matrix:
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
    %   2) Y - N-by-1 vector of solution at time t
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


























