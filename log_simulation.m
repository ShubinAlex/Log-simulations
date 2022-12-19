function [simulated_log] = log_simulation( log, z, sigma, L, nugget, range, N)
% LOG_SIMULATION variogram-based log simulations
%   log - log curve values (e.g. porosty, GR ...)
%   z - depth values
%   sigma - smoothing parameter for trend extraction (10-20)
%   L - max number of lags (lag - distance between neighbor points)
%   nugget - nugget (%) value for spherical variogram
%   range - range (m) value for spherical variogram
%   N - number of simulations
%   OUTPUT
%   simulated logs matrix (z*N)

%   Alexey Shubin 2022. Reference Dvorkin J. Seismic reflection of rock
%   properties, 2014 p. 95
%   see lusimulation.m


simulated_log = zeros(length(z),N);


% add samples to log data
up_matrix = zeros(500,1);
up_matrix(1:500)=log(1);
down_matrix = zeros(500,1);
down_matrix(1:500) = log(max(size(log)));
extended_log = [up_matrix;log];
extended_log = [extended_log;down_matrix];


log_filtered = gaussfilt(z,extended_log,sigma); % filtering

log_filtered = log_filtered(501:500+max(size(log)));
x = log-log_filtered;

for j = 1:N
    
    w = lusimulation(z, x, L, nugget, range); % Cholesky decomposition
    
    simulated_log(:,j)=log_filtered+w;
    simulated_log(simulated_log<0)=0;
    subplot(2,1,2); plot(z, log, 'k', 'linewidth', 2)
    grid on;
    hold on;
    plot(z, simulated_log(:,j), 'color', [.5 .5 .5]);
    plot(z, log_filtered, 'g');
    legend('initial log','simulated log', 'trend');
    xlabel('depth, m');
    ylabel('log');
    hold on
end

end

