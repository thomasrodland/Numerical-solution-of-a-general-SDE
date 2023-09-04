function Z = normal_generator(N, m1, m2)
    % N: Number of random samples to generate
    % m1, m2: Mean and standard deviation of the normal distribution

    % Initialize an array to store the generated random numbers
    Z = zeros(1, N);

    % Generate N/2 pairs of independent standard normal random numbers
    for i = 1:N/2
        % Generate two independent uniform random numbers in (0,1)
        u1 = rand;
        u2 = rand;

        % Use the Box-Muller transform to generate two independent
        % standard normal random numbers
        z1 = sqrt(-2 * log(u1)) * cos(2 * pi * u2);
        z2 = sqrt(-2 * log(u1)) * sin(2 * pi * u2);

        % Apply the mean and standard deviation transformations
        Z(2*i-1) = m1 + z1;
        Z(2*i) = m2 + z2;
    end

    % If N is odd, generate one more random number and append it
    if mod(N, 2) == 1
        u = rand;
        z = sqrt(-2 * log(u)) * cos(2 * pi * rand);
        Z(end+1) = m1 + z;
    end
end
