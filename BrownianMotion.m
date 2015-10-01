function B = BrownianMotion(dt,K,M)

    if nargin < 3
        M = 1;
    end
    
    sqrt_dt = sqrt(dt);

	B = zeros(M,K);
    for i = 2:K
        B(:,i) = B(:,i-1) + sqrt_dt * randn(M,1);
    end

end