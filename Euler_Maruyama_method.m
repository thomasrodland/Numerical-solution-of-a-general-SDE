function [lt, X_EM] = Euler_Maruyama_method(a, b, T, N, x0, m1, m2)
    % Parameters:
    % a: Drift coefficient function a(t, x)
    % b: Diffusion coefficient function b(t, x)
    % T: Total simulation time
    % N: Number of time steps
    % x0: Initial value
    % m1, m2: Mean values for normal_generator
    
    % Initialize arrays to store time and solution
    lt = linspace(0, T, N + 1);
    dt = T / N;
    X_EM = zeros(1, N + 1);
    X_EM(1) = x0;

    % Generate N independent standard normal random numbers
    Z = normal_generator(N, m1, m2);

    % Perform Euler-Maruyama simulation
    for i = 1:N
        dW = Z(i);
        X_EM(i + 1) = X_EM(i) + a(lt(i), X_EM(i)) * dt + b(lt(i), X_EM(i)) * dW;
    end
end
