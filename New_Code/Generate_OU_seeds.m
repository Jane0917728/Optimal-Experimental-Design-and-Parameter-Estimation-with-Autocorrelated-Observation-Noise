function OUrandom=Generate_OU_seeds(phi,Sigma,t,Z)


    % phi: Mean reversion rate
    % Sigma: Diffusion coefficient
    % t: Time vector
    % Z: Sequence of standard normal random variables, length n_s
    n_s = length(t); % Number of time points

    % Initialize the OU process sequence
    X = zeros(1, n_s); 

    % Initialize X(1) based on randomness, X(1) ~ N(0, sigma^2)
    X(1) = Sigma * Z(1);  % The first element is random

    % Compute the values of the OU process at each time point
    for i = 1:n_s-1
        dt = t(i+1) - t(i); % Compute the time step
        mean_X = X(i) * exp(-phi * dt);  % Mean reversion term
        variance_X = (Sigma^2 / (2 * phi)) * (1 - exp(-2 * phi * dt));  % Variance term
        
        % Generate increment from normal distribution
        X(i+1) = mean_X + sqrt(variance_X) * Z(i+1);
    end
    OUrandom=X;
    % figure;
    %   plot(t_opt, X, '+m-',t_opt, X1, 'dg','MarkerSize',6,'LineWidth', 1.2);
end


 
