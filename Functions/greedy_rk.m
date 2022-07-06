function [x_resp, sampled_freq] = greedy_rk(M, K, F, freq, tol)

x_resp = [];
sampled_freq = [freq(1)];
%% initialize the reduced system
omega = 2*pi*freq(1); 
V = compute_basis(M, K, F, omega);
[M_red, K_red, F_red] = compute_projection(M, K, F, V);

%% Frequency sweep
for iFreq = 1:numel(freq)
    omega = 2*pi*freq(iFreq);
    x = compute_response(M_red, K_red, F_red, V, omega);
    error = norm((-omega^2*M + K)*x-F)
    if error > tol
        sampled_freq = [sampled_freq, freq(iFreq)]
        V = [V, compute_basis(M, K, F, omega)];
        V = full(V); 
        V = orth(V); 
        [M_red, K_red, F_red] = compute_projection(M, K, F, V);
        x = compute_response(M_red, K_red, F_red, V, omega);
    end
    x_resp = [x_resp, x];
end

function v = compute_basis(M, K, F, omega)
G = -omega^2*M + K;
dG = decomposition(G,'lu');
v = dG\F;

function [M_red, K_red, F_red] = compute_projection(M, K, F, V)
M_red = V'*M*V;
K_red = V'*K*V;
F_red = V'*F;

function x = compute_response(M_red, K_red, F_red, V, omega)
G_red = -omega^2*M_red + K_red;
dG_red = decomposition(G_red,'lu');
x = V*(dG_red\F_red);

function error_freq = compute_error(M_red, K_red, F_red, V, M, K, F, freq)
error_freq = [];
freq = freq(1:10:end);
for iFreq = 1:numel(freq)
    iFreq
    omega = 2*pi*freq(iFreq);
    x = V*((-omega^2*M_red + K_red)\F_red);
    error_freq(end+1) = norm((-omega^2*M + K)*x-F);
end
