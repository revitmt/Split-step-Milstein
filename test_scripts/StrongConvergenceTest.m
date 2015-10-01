function [h,err,order] = StrongConvergenceTest(method,DriftVector,DiffusionMatrix,FTime,Y0)

%   Purpose
%   =======
%   Test strong convergence of the numerical method
%
%   Method
%   ======
%   Strong error is obtained using the reference solution calculated
%   with the Milstein method on the fine reference grid
%
%   IN
%   ==
%   1) method          - function handle that calls the corresponding
%                        numerical method
%   2) DriftVector     - function handle that evaluates drift vector
%   3) DiffusionMatrix - function handle that evaluates matrix of 
%                        diffusion coefficients
%   4) FTime - final time of the simulation
%   5) Y0 - n-dimensional column vector with initial data
%
%
%   OUT
%   ===
%   h     - vector with time steps
%   err   - vector with corresponding strong errors
%   order - estimated order of convergence


    subdiv    = 15;
    dt_ref    = 2^(-subdiv);
    tspan_ref = 0:dt_ref:FTime;   % reference time discretization
    tspan_ref(end) = tspan_ref(end-1) + dt_ref;

    M     = size(DiffusionMatrix(1,Y0),2);  % number of noise channels
    K_ref = length(tspan_ref);              % number of time points
    Nm    = length(method);                 % number of methods to test

    test_lev_1 = 1;               % initial discretization
    test_lev_2 = subdiv-10;        % final discretization

    NPaths = 50;                  % number of simulated paths

    h   = zeros( 1,  test_lev_2 - test_lev_1 + 1 );
    err = zeros( Nm, test_lev_2 - test_lev_1 + 1 );

    for path = 1:NPaths  
        % reference solution on the fine grid
        Wiener_ref = BrownianMotion(dt_ref,K_ref,M);
        Y_ref = Milstein(DriftVector,DiffusionMatrix,tspan_ref,Y0,Wiener_ref);

        % run test levels
        for p = test_lev_1:test_lev_2
            mult = 2^p;
            dt = mult * dt_ref;
            tspan = 0:dt:FTime;
            K = length(tspan);
            Wiener = zeros(M,K);
            for i = 1:K
                j = (i-1)*mult + 1;
                Wiener(:,i) = Wiener_ref(:,j);
            end
    
            % test different methods
            for n = 1:Nm
                Y = method{n}(DriftVector,DiffusionMatrix,tspan,Y0,Wiener);
                err(n, p-test_lev_1+1) = err(n, p-test_lev_1+1) + norm(Y_ref(:,end)-Y(:,end));
            end
            h(p-test_lev_1+1) = dt;
    
            disp(sprintf('Path: %i; Level: %i',path,p));
        end

    end

    err = err ./ NPaths;

    order = zeros(Nm, 1);
    for n = 1:Nm
        p = polyfit( log(h), log(err(n,:)), 1 );
        order(n) = p(1);
    end

end