clear;
clc;

T=1;
x0 = 10;
mu = 0.5;
sigma = 0.3;
a = @(t,x) mu * x;
b = @(t,x) sigma * x;
diff_b = @(t,x) sigma;
m1 = 0.2521;
m2 = 0.14324;

h_values = 0.005 * (1 / 2).^ (0:3);
num_simulations = 500000;

orders_strong_EM = zeros(size(h_values));
orders_weak_EM = zeros(size(h_values));
orders_strong_Milstein = zeros(size(h_values));
orders_weak_Milstein = zeros(size(h_values));

errors_Milstein_weak = zeros(size(h_values));
errors_Milstein_strong = zeros(size(h_values));
errors_EM_weak = zeros(size(h_values));
errors_EM_strong = zeros(size(h_values));

for i = 1:length(h_values)
    h = h_values(i);
    N = T / h;
    errors_EM_w = zeros(1, num_simulations);
    errors_Milstein_w = zeros(1, num_simulations);
    errors_EM_s = zeros(1, num_simulations);
    errors_Milstein_s = zeros(1, num_simulations);
    history_S = zeros(1, num_simulations);

    for j = 1:num_simulations
        m1 = m1 - m1/num_simulations;
        m2 = m2 - m2/num_simulations;

        Z=normal_generator(N,m1,m2);
        [lt,X_EM] = Euler_Maruyama_method(a,b,T,N,x0,m1,m2);
        [lt,X_M] = Milstein_method(a,b,diff_b,T,N,x0,m1,m2);

        Z=normal_generator(N,m1,m2);
        S = zeros(size(X_M));
        S(1)=x0;
        B = [0; cumsum(Z')] *sqrt(h);
        for k=1:N
          S(k+1) = x0 * exp((mu - (sigma^2)/2) * lt(k+1) + sigma*B(k+1));
        end

        errors_EM_w(j) = X_EM(end);
        errors_Milstein_w(j) = X_M(end);
        history_S(j) = S(end);
    end
    disp(i)
    errors_Milstein_weak(i) = abs(mean(errors_Milstein_w - history_S));
    errors_EM_weak(i) = abs(mean(errors_EM_w - history_S));
    errors_Milstein_strong(i) = mean(abs(errors_Milstein_w - history_S));
    errors_EM_strong(i) = mean(abs(errors_EM_w - history_S));
end

errors_log_Milstein_strong = log(errors_Milstein_strong);
errors_log_EM_strong = log(errors_EM_strong);
errors_log_Milstein_weak = log(errors_Milstein_weak);
errors_log_EM_weak = log(errors_EM_weak);

log_h_values = log(h_values);
coefficients_Milstein_strong = polyfit(log_h_values, errors_log_Milstein_strong, 1)
coefficients_EM_strong = polyfit(log_h_values, errors_log_EM_strong, 1)
coefficients_Milstein_weak = polyfit(log_h_values, errors_log_Milstein_weak, 1)
coefficients_EM_weak = polyfit(log_h_values, errors_log_EM_weak, 1)

x = linspace(min(log_h_values), max(log_h_values), 100);

fit_Milstein_strong = polyval(coefficients_Milstein_strong, x);
fit_EM_strong = polyval(coefficients_EM_strong, x);

fit_Milstein_weak = polyval(coefficients_Milstein_weak, x);
fit_EM_weak = polyval(coefficients_EM_weak, x);

%Create separate figures as tiles for each fit
figure;
tiledlayout(2, 2);

nexttile;
scatter(log_h_values, errors_log_Milstein_strong, 'DisplayName', 'Milstein Strong');
hold on;
plot(x, fit_Milstein_strong, 'DisplayName', 'Milstein Strong Fit');
xlabel('log(h) values');
ylabel('Errors');
title('Milstein Strong Fit');
legend('Location', 'best');

nexttile;
scatter(log_h_values, errors_log_EM_strong, 'DisplayName', 'EM Strong');
hold on;
plot(x, fit_EM_strong, 'DisplayName', 'EM Strong Fit');
xlabel('log(h) values');
ylabel('Errors');
title('EM Strong Fit');
legend('Location', 'best');

nexttile;
scatter(log_h_values, errors_log_Milstein_weak, 'DisplayName', 'Milstein Weak');
hold on;
plot(x, fit_Milstein_weak, 'DisplayName', 'Milstein Weak Fit');
xlabel('log(h) values');
ylabel('Errors');
title('Milstein Weak Fit');
legend('Location', 'best');

nexttile;
scatter(log_h_values, errors_log_EM_weak, 'DisplayName', 'EM Weak');
hold on;
plot(x, fit_EM_weak, 'DisplayName', 'EM Weak Fit');
xlabel('log(h) values');
ylabel('Errors');
title('EM Weak Fit');
legend('Location', 'best');

%Adjust the layout and spacing of the tiles
padding = 0.03;
margin = 0.05;
tiledlayout.Padding = 'compact';
tiledlayout.TileSpacing = 'compact';
tiledlayout.Padding = padding;
tiledlayout.Margin = margin;