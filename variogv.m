function [ gamma ] = variogv( z, x, L, nugget, range )
%variogv vertical (1D) variogram calculation and model estimation
%   z - spatial parameter (e.g. depth, m)
%   x - variable (e.g. porosity)
%   L - max number of lags (lag - distance between neighbor points)
%   model variogram parameters:
%   nugget - in percents(%) from sill
%   range - in meters (m)
%   OUTPUT
%   gamma  - variogram values


if nargin==3,
nugget = 0;
range = 2;
end

% lag value
delta_z = z(2)-z(1);
% lag vector
lag_vector = 1:L;
lag_distance = delta_z*lag_vector;
% data length
N = length(z);
% experimental variogram
count_pairs = zeros(1,L);
gamma = zeros(1,L);
% calculation (x(j+i)-x(j))^2 for pairs with lags from 1 to L
pairs = zeros(N-1,L);

for i = 1:L
    for j = 1:N-i
        pairs(j,i) = (x(j+i)-x(j))^2;
        count_pairs(i) = N-i;
    end

gamma(i) = 0.5*(sum(pairs(:,i)/count_pairs(i)));

end

scatter(lag_distance,gamma); % plot experimental variogram
grid on;
xlabel('distance, m');
ylabel('gamma');
hold on

% variogram model
sill = var(x(1:L));
nugget = (nugget/100)*sill;
sph_var=zeros(1,L); % spherical variogram
for i=1:L
    if lag_distance(1,i) <= range
        sph_var(1,i)= nugget+sill*(1.5*(lag_distance(1,i)/range)-0.5*(lag_distance(1,i)/range)^3);
    else
        sph_var(1,i)= nugget+sill;
    end
end
plot(lag_distance, sph_var, 'r');


end

