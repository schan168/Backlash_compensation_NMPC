clear; close all; clc;

params.Jm = 0.02;
params.Jl = 0.05;
params.Bm = 0.01;
params.Bl = 0.01;
params.Ks = 5;
params.Bs = 0.05;
params.d  = 0.03;

Ts = 0.01;
Np = 50;
Nc = 8;

umin = -1;
umax = 1;

Q = diag([5, 0.5, 300, 5]);
R = 0.1;
Rd = 5;
Q_terminal = 3;

Tsim = 6;
Nsim = round(Tsim/Ts);
t = (0:Nsim)'*Ts;

theta_l_ref = 0.5 * ones(Nsim+1, 1);

x_hist = zeros(4, Nsim+1);
u_hist = zeros(Nsim, 1);
x_hist(:,1) = [0; 0; 0; 0];

u_prev = 0;
u_warm = zeros(Nc, 1);

n_sub = 10;
dt = Ts / n_sub;

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
                       'MaxIterations', 150, 'OptimalityTolerance', 1e-5, ...
                       'StepTolerance', 1e-8);

fprintf('Running Optimized NMPC...\n');

tic;
for k = 1:Nsim
    x_current = x_hist(:,k);
    
    ref_horizon = zeros(4, Np);
    for j = 1:Np
        idx = min(k+j, Nsim+1);
        target_load = theta_l_ref(idx);
        target_motor = target_load + params.d + 0.005;
        ref_horizon(:,j) = [target_motor; 0; target_load; 0];
    end
    
    costFun = @(u_seq) nmpcCost(u_seq, x_current, u_prev, ref_horizon, params, Ts, Np, Nc, Q, R, Rd, Q_terminal);
    
    lb = umin * ones(Nc, 1);
    ub = umax * ones(Nc, 1);
    
    A_rate = zeros(Nc-1, Nc);
    for i = 1:Nc-1
        A_rate(i,i) = -1;
        A_rate(i,i+1) = 1;
    end
    b_rate_upper = 0.05 * ones(Nc-1, 1);
    A_ineq = [A_rate; -A_rate];
    b_ineq = [b_rate_upper; b_rate_upper];
    
    try
        u_opt = fmincon(costFun, u_warm, A_ineq, b_ineq, [], [], lb, ub, [], options);
        if isempty(u_opt)
            u_opt = u_warm;
        end
    catch
        u_opt = u_warm;
    end
    
    u = u_opt(1);
    u_hist(k) = u;
    u_prev = u;
    u_warm = [u_opt(2:end); u_opt(end)];
    
    x_temp = x_current;
    for s = 1:n_sub
        x_temp = nlDynamics(x_temp, u, dt, params);
    end
    x_hist(:,k+1) = x_temp;
    
    if mod(k, 100) == 0
        fprintf('Step %d/%d | theta_l = %.4f | u = %.4f\n', k, Nsim, x_hist(3,k+1), u);
    end
end
comp_time = toc;

fprintf('Completed in %.2f seconds (%.2f ms/step)\n\n', comp_time, comp_time/Nsim*1000);

d = params.d;
gap = x_hist(1,:) - x_hist(3,:);
in_gap = abs(gap) <= d;
load_error = x_hist(3,:)' - theta_l_ref;

rms_error = sqrt(mean(load_error.^2));
ss_error = abs(mean(x_hist(3,end-50:end)) - theta_l_ref(end));
time_engaged = sum(~in_gap)/length(in_gap)*100;

final_ref = theta_l_ref(end);
settled_idx = find(abs(x_hist(3,:) - final_ref) < 0.02*abs(final_ref), 1, 'first');
if isempty(settled_idx)
    settling_time = Inf;
else
    settling_time = settled_idx * Ts;
end

overshoot = max(0, (max(x_hist(3,:)) - final_ref)/final_ref*100);

fprintf('==================== RESULTS ====================\n');
fprintf('RMS Load Error:      %.4f rad\n', rms_error);
fprintf('Steady-State Error:  %.4f rad\n', ss_error);
fprintf('Settling Time:       %.3f s\n', settling_time);
fprintf('Overshoot:           %.2f %%\n', overshoot);
fprintf('Time Engaged:        %.1f %%\n', time_engaged);
fprintf('Final Position:      %.4f rad (ref: %.4f)\n', x_hist(3,end), final_ref);
fprintf('=================================================\n');

figure('Name', 'Optimized NMPC Results', 'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.9]);

subplot(4,1,1);
plot(t, x_hist(3,:), 'b-', 'LineWidth', 2); hold on;
plot(t, theta_l_ref, 'r--', 'LineWidth', 1.5);
plot(t, x_hist(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
grid on;
ylabel('\theta (rad)');
legend('\theta_l (Load)', 'Reference', '\theta_m (Motor)', 'Location', 'southeast');
title(sprintf('Position Tracking | SS Error: %.4f rad | Settling: %.3fs', ss_error, settling_time));

subplot(4,1,2);
stairs(t(1:end-1), u_hist, 'b-', 'LineWidth', 1.5); hold on;
yline(umin, 'r--'); yline(umax, 'r--');
grid on;
ylabel('u (N·m)');
title('Control Input');
ylim([umin-0.2 umax+0.2]);

subplot(4,1,3);
plot(t, gap, 'b-', 'LineWidth', 1.5); hold on;
yline(d, 'r--', 'LineWidth', 1.5);
yline(-d, 'r--', 'LineWidth', 1.5);
fill([t(1) t(end) t(end) t(1)], [-d -d d d], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on;
ylabel('Gap (rad)');
title(sprintf('Backlash Gap (Dead-zone: ±%.3f rad) | Engaged: %.1f%%', d, time_engaged));

subplot(4,1,4);
plot(t, load_error, 'b-', 'LineWidth', 1.5); hold on;
yline(0, 'k--');
fill([t(1) t(end) t(end) t(1)], [-0.02*final_ref -0.02*final_ref 0.02*final_ref 0.02*final_ref], ...
     'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
grid on;
ylabel('Error (rad)');
xlabel('Time (s)');
title(sprintf('Tracking Error | RMS: %.4f rad | Overshoot: %.2f%%', rms_error, overshoot));

sgtitle('Nonlinear MPC - Two Mass System with Backlash', 'FontSize', 14, 'FontWeight', 'bold');

function J = nmpcCost(u_seq, x0, u_prev, ref_horizon, params, Ts, Np, Nc, Q, R, Rd, Q_terminal)
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
            x = nlDynamics(x, u, dt, params);
        end
        
        e = x - ref_horizon(:,j);
        
        if j < Np
            J = J + e' * Q * e;
        else
            J = J + Q_terminal * (e' * Q * e);
        end
        
        if j <= Nc
            J = J + R * u^2 + Rd * (u - u_last)^2;
            u_last = u;
        end
    end
end

function x_new = nlDynamics(x, u, dt, params)
    Jm = params.Jm;
    Jl = params.Jl;
    Bm = params.Bm;
    Bl = params.Bl;
    Ks = params.Ks;
    Bs = params.Bs;
    d = params.d;
    
    theta_m = x(1);
    omega_m = x(2);
    theta_l = x(3);
    omega_l = x(4);
    
    rel = theta_m - theta_l;
    rel_dot = omega_m - omega_l;
    
    if abs(rel) <= d
        tau_s = 0;
    else
        tau_s = Ks * (rel - sign(rel)*d) + Bs * rel_dot;
    end
    
    alpha_m = (u - Bm*omega_m - tau_s) / Jm;
    alpha_l = (tau_s - Bl*omega_l) / Jl;
    
    x_new = zeros(4,1);
    x_new(1) = theta_m + dt*omega_m + 0.5*dt^2*alpha_m;
    x_new(2) = omega_m + dt*alpha_m;
    x_new(3) = theta_l + dt*omega_l + 0.5*dt^2*alpha_l;
    x_new(4) = omega_l + dt*alpha_l;
end
