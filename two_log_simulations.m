function [sim_log1, sim_log2 ] = two_log_simulations(log1, log2, z, sigma, nugget, range1, range2, N)
% TWO_LOG_SIMULATIONS - stochastic simulation of related variables
%   log1 - first log curve values (e.g. porosty, GR ...)
%   log2 - second log curve values (e.g. porosty, GR ...)
%   z - depth values
%   sigma - smoothing parameter for trend extraction (gaussian smoothing) (10-20)
%   
%   nugget - nugget value (%) for spherical variogram
%   range - range value (m) for spherical variogram
%   N - number of simulations
% OUTPUT
%   simulated logs matrix (z*N)

%   Alexey Shubin 2022. Reference Dvorkin J. Seismic reflection of rock
%   properties, 2014 p. 98

L = length(z); %max number of lags (lag - distance between neighbor points)
sim_log1 = zeros(length(z),N);
sim_log2 = zeros(length(z),N);

% add samples to log1 data
up_matrix=zeros(500,1);
up_matrix(1:500)=log1(1);
down_matrix=zeros(500,1);
down_matrix(1:500)=log1(max(size(log1)));
extended_log1=[up_matrix;log1];
extended_log1=[extended_log1;down_matrix];

log_filtered1=gaussfilt(z,extended_log1,sigma);
log_filtered1=log_filtered1(501:500+max(size(log1)));
x=log1-log_filtered1;

% add samples to log2 data
up_matrix=zeros(500,1);
up_matrix(1:500)=log2(1);
down_matrix=zeros(500,1);
down_matrix(1:500)=log2(max(size(log2)));
extended_log2=[up_matrix;log2];
extended_log2=[extended_log2;down_matrix];

log_filtered2=gaussfilt(z,extended_log2,sigma);
log_filtered2=log_filtered2(501:500+max(size(log2)));
y=log2-log_filtered2;

subplot(1,2,1)
[~, sph_var1, ~]=variogv(z, x, L, nugget, range1);
sph_var1=sph_var1./var(x);

subplot(1,2,2)
[~, sph_var2, ~]=variogv(z, y, L, nugget, range2);
sph_var2=sph_var2./var(y);

hold on;

sph_var=0.5.*(sph_var1+sph_var2);

Cov_model=1-sph_var; % covariance model calculation
J=max(size(log1));
Z=zeros(1,J-L); % add 0's to log size
Cov_model=[Cov_model, Z];
Cov_model_mat=toeplitz(Cov_model); % covariance model matrix

S=cov(x,y,1); % data covariance

Kron_prod=kron(S,Cov_model_mat);

R=chol(Kron_prod);

figure;
hold on;
subplot(1,2,1)

plot(log1, z, 'k', 'linewidth', 2);
set(gca, 'YDir','reverse')
grid on;
hold on;

subplot(1,2,2)
plot(log2, z, 'k', 'linewidth', 2);
set(gca, 'YDir','reverse')
grid on;
hold on;

for j = 1:N

u=random('Normal',0,1,2.*J,1); % u - uncorrelated vector
w=R*u;
sim=[log_filtered1;log_filtered2]+w;
sim_log1(:,j)=sim(1:J);
sim_log2(:,j)=sim(J+1:2*J);

sim_log1(sim_log1<0)=0;
sim_log2(sim_log2<0)=0;


subplot(1,2,1)

plot(sim_log1(:,j), z, 'color', [.5 .5 .5]);
xlabel('log1');
ylabel('depth, m');

subplot(1,2,2)

plot(sim_log2(:,j), z, 'color', [.5 .5 .5]);
xlabel('log2');
ylabel('depth, m');
hold on;

end

end

