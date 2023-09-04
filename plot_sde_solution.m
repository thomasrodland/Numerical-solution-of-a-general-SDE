
% Define your SDE parameters and functions
a = @(t, x) mu * x;
b = @(t, x) sigma * x;
diff_b = @(t, x) sigma;
T = 1;
N = 1000;
x0 = 10;
m1 = 0.2521;
m2 = 0.14324;

% Plot the SDE solution using the Euler method
plot_sde_solution2('Euler', a, b, diff_b, T, N, x0, m1, m2);

function plot_sde_solution2(method, a, b, diff_b, T, N, x0, m1, m2)
    % Parameters:
    % method: 'Euler' or 'Milstein' to specify the simulation method
    % a: Drift coefficient function a(t, x)
    % b: Diffusion coefficient function b(t, x)
    % diff_b: Derivative of diffusion coefficient function b(t, x)
    % T: Total simulation time
    % N: Number of time steps
    % x0: Initial value
    % m1, m2: Mean values for normal_generator
    
    % Simulate the SDE using the specified method
    if strcmp(method, 'Euler')
        [lt, X] = Euler_Maruyama_method(a, b, T, N, x0, m1, m2);
    elseif strcmp(method, 'Milstein')
        [lt, X] = Milstein_method(a, b, diff_b, T, N, x0, m1, m2);
    else
        error('Invalid simulation method. Use ''Euler'' or ''Milstein''.');
    end

    % Plot the solution
    figure;
    plot(lt, X);
    xlabel('Time');
    ylabel('Solution');
    title(['SDE Solution using ' method ' Method']);
    grid on;
end


