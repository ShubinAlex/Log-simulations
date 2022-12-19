function w=lusimulation(z, x, L, nugget, range)
% lusimulation   simulate a random vertically correlated vector w
    
%   z - spatial parameter (e.g. depth, m)
%   x - variable (e.g. porosity)
%   L - max number of lags (lag - distance between neighbor points)
%   model variogram parameters:
%   nugget - in percents(%) from sill
%   range - in meters (m)
%   OUTPUT
%   w  - random vertically correlated vector

if nargin==3,
nugget = 0;
range = 2;
end

% lag value
delta_z = z(2)-z(1);
% lag vector
lag_vector = 1:L;
lag_distance = delta_z*lag_vector;
S=max(size(z));
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

subplot(2,1,1);
scatter(lag_distance,gamma, 'k'); % plot experimental variogram
grid on;
title('variogram');
xlabel('distance, m');
ylabel('gamma');
hold on

% variogram model
sill=var(x(1:L));
nugget = (nugget/100)*sill;
sph_var=zeros(1,L); % spherical variogram
for i=1:L
    if lag_distance(1,i) <= range
        sph_var(1,i)=nugget+sill*(1.5*(lag_distance(1,i)/range)-0.5*(lag_distance(1,i)/range)^3);
    else
        sph_var(1,i)= nugget+sill;
    end
end

plot(lag_distance, sph_var, 'r');
hold on
%   simulate a random vertically correlated vector w
   
   Cov_model=var(x(1:L))-sph_var; % covariance model calculation
   Cov_model(Cov_model<0)=0; % remove negative values
   Z=zeros(1,S-L); % add 0's to size z
   Cov_model=[Cov_model,Z];
   Cov_model_mat=toeplitz(Cov_model); % covariance matrix
   R=chol(Cov_model_mat, 'lower'); % Cholesky decomposition
   u=random('Normal',0,1,S,1); % u - uncorrelated vector
   w=R*u; % simulated vector
   
end