function [lt, X_M] = Milstein_method(a, b, diff_b, T, N, x0, m1, m2)
    % Parameters:
    % a: Drift coefficient function a(t, x)
    % b: Diffusion coefficient function b(t, x)
    % diff_b: Derivative of diffusion coefficient function b(t, x)
    % T: Total simulation time
    % N: Number of time steps
    % x0: Initial value
    % m1, m2: Mean values for normal_generator
    
    % Initialize arrays to store time and solution
    lt = linspace(0, T, N + 1);
    dt = T / N;
    X_M = zeros(1, N + 1);
    X_M(1) = x0;

    % Generate N independent standard normal random numbers
    Z = normal_generator(N, m1, m2);

    % Perform Milstein simulation
    for i = 1:N
        dW = Z(i);
        X_M(i + 1) = X_M(i) + a(lt(i), X_M(i)) * dt + b(lt(i), X_M(i)) * dW + ...
                     0.5 * b(lt(i), X_M(i)) * diff_b(lt(i), X_M(i)) * (dW^2 - dt);
    end
end
