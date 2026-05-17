clear; close all; clc;

params.Jm = 0.02;
params.Jl = 0.05;
params.Bm = 0.01;
params.Bl = 0.01;
params.Ks = 5;
params.Bs = 0.05;
params.d  = 0.03;

Ts = 0.01;
Np = 30;
Nc = 10;

umin = -1;
umax = 1;
dumax = 0.1;

Q = diag([10, 0.1, 200, 1]);
R = 0.001;
Rd = 0.01;
Q_dz = 50;

Tsim = 6;
Nsim = round(Tsim/Ts);
t = (0:Nsim)'*Ts;

theta_l_ref = 0.5 * ones(Nsim+1, 1);

x_hist = zeros(4, Nsim+1);
u_hist = zeros(Nsim, 1);
x_hist(:,1) = [0; 0; 0; 0];

u_prev = 0;
u_warm = 0.5 * ones(Nc, 1);

n_sub = 10;
dt = Ts / n_sub;

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
                       'MaxIterations', 100, 'OptimalityTolerance', 1e-4);

fprintf('Running NMPC \n');

tic;
for k = 1:Nsim
    x_current = x_hist(:,k);
    
    ref_horizon = zeros(4, Np);
    for j = 1:Np
        idx = min(k+j, Nsim+1);
        ref_horizon(:,j) = [theta_l_ref(idx) + params.d*1.5; 0; theta_l_ref(idx); 0];
    end
    
    costFun = @(u_seq) nmpcCost_DeadzonePenalty(u_seq, x_current, u_prev, ref_horizon, params, Ts, Np, Nc, Q, R, Rd, Q_dz);
    
    lb = umin * ones(Nc, 1);
    ub = umax * ones(Nc, 1);
    
    try
        u_opt = fmincon(costFun, u_warm, [], [], [], [], lb, ub, [], options);
        if isempty(u_opt)
            u_opt = u_warm;
        end
    catch
        u_opt = u_warm;
    end
    
    u = max(min(u_opt(1), umax), umin);
    u_hist(k) = u;
    u_prev = u;
    u_warm = [u_opt(2:end); u_opt(end)];
    
    x_temp = x_current;
    for s = 1:n_sub
        x_temp = nonlinearDynamics(x_temp, u, dt, params);
    end
    x_hist(:,k+1) = x_temp;
    
    if mod(k, 100) == 0
        fprintf('Step %d/%d | theta_l = %.4f | u = %.4f\n', k, Nsim, x_hist(3,k+1), u);
    end
end
comp_time = toc;

fprintf('Completed in %.2f seconds (%.2f ms/step)\n', comp_time, comp_time/Nsim*1000);

d = params.d;
gap = x_hist(1,:) - x_hist(3,:);
in_gap = abs(gap) <= d;
load_error = x_hist(3,:)' - theta_l_ref;

rms_error = sqrt(mean(load_error.^2));
ss_error = mean(load_error(end-100:end));
time_engaged = sum(~in_gap)/length(in_gap)*100;

final_ref = theta_l_ref(end);
settled = abs(x_hist(3,:) - final_ref) < 0.02*abs(final_ref);
idx = find(settled, 1, 'first');
if isempty(idx)
    settling_time = Inf;
else
    settling_time = idx * Ts;
end

overshoot = max(0, (max(x_hist(3,:)) - final_ref)/final_ref*100);
IAE = sum(abs(load_error))*Ts;

fprintf('\nResults\n');
fprintf('  RMS Load Error:      %.4f rad\n', rms_error);
fprintf('  Steady-State Error:  %.4f rad\n', ss_error);
fprintf('  Time Engaged:        %.1f %%\n', time_engaged);
fprintf('  Settling Time:       %.3f s\n', settling_time);
fprintf('  Overshoot:           %.2f %%\n', overshoot);
fprintf('  Final Load Position: %.4f rad (ref: %.4f)\n', x_hist(3,end), final_ref);

figure('Name', 'NMPC with Dead-Zone Penalty', 'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.9]);

subplot(5,1,1);
plot(t, x_hist(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(t, x_hist(3,:), 'r-', 'LineWidth', 1.5);
plot(t, theta_l_ref, 'k--', 'LineWidth', 1.5);
grid on; ylabel('\theta (rad)');
legend('\theta_m', '\theta_l', 'Reference', 'Location', 'southeast');
title('Positions');

subplot(5,1,2);
stairs(t(1:end-1), u_hist, 'b-', 'LineWidth', 1.5); hold on;
yline(umin, 'r--'); yline(umax, 'r--');
grid on; ylabel('u (N·m)');
title('Control Input');

subplot(5,1,3);
plot(t, gap, 'b-', 'LineWidth', 1.5); hold on;
yline(d, 'r--', 'LineWidth', 1.5);
yline(-d, 'r--', 'LineWidth', 1.5);
patch([t(1) t(end) t(end) t(1)], [-d -d d d], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
grid on; ylabel('Gap (rad)');
title('Relative Displacement');

subplot(5,1,4);
area(t, double(~in_gap), 'FaceColor', [0.3 0.7 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
area(t, double(in_gap), 'FaceColor', [0.9 0.3 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
ylim([-0.1 1.1]);
yticks([0 1]); yticklabels({'Dead-Zone', 'Engaged'});
grid on;
title(sprintf('Engagement Status (%.1f%% engaged)', time_engaged));

subplot(5,1,5);
plot(t, load_error, 'b-', 'LineWidth', 1.5);
grid on; ylabel('Error (rad)'); xlabel('Time (s)');
title(sprintf('Tracking Error (RMS: %.4f rad)', rms_error));

sgtitle('NMPC with Dead-Zone Penalty', 'FontSize', 14, 'FontWeight', 'bold');

function J = nmpcCost_DeadzonePenalty(u_seq, x0, u_prev, ref_horizon, params, Ts, Np, Nc, Q, R, Rd, Q_dz)
    d = params.d;
    n_sub = 5;
    dt = Ts/n_sub;
    
    J = 0;
    x = x0;
    u_last = u_prev;
    
    for j = 1:Np
        if j <= Nc
            u = u_seq(j);
        else
            u = u_seq(Nc);
        end
        
        for s = 1:n_sub
            x = nonlinearDynamics(x, u, dt, params);
        end
        
        e = x - ref_horizon(:,j);
        J = J + e' * Q * e;
        
        gap = x(1) - x(3);
        if abs(gap) < d
            J = J + Q_dz * (d - abs(gap))^2;
        end
        
        if j <= Nc
            J = J + R * u^2 + Rd * (u - u_last)^2;
            u_last = u;
        end
    end
end

function x_new = nonlinearDynamics(x, u, dt, params)
    Jm = params.Jm; Jl = params.Jl;
    Bm = params.Bm; Bl = params.Bl;
    Ks = params.Ks; Bs = params.Bs;
    d = params.d;
    
    theta_m = x(1); theta_m_dot = x(2);
    theta_l = x(3); theta_l_dot = x(4);
    
    rel = theta_m - theta_l;
    rel_dot = theta_m_dot - theta_l_dot;
    
    if abs(rel) <= d
        tau_s = 0;
    else
        tau_s = Ks * (rel - sign(rel)*d) + Bs * rel_dot;
    end
    
    theta_m_ddot = (u - Bm*theta_m_dot - tau_s) / Jm;
    theta_l_ddot = (tau_s - Bl*theta_l_dot) / Jl;
    
    x_new = [theta_m + dt*theta_m_dot;
             theta_m_dot + dt*theta_m_ddot;
             theta_l + dt*theta_l_dot;
             theta_l_dot + dt*theta_l_ddot];
end
