clear; close all; clc;

params.Jm = 0.02;
params.Jl = 0.05;
params.Bm = 0.01;
params.Bl = 0.01;
params.Ks = 5;
params.Bs = 0.05;
params.d  = 0.03;

Ts = 0.01;
Tsim = 6;
Nsim = round(Tsim/Ts);
t = (0:Nsim)'*Ts;

umin = -1;
umax = 1;
dumax = 0.1;

theta_l_ref = 0.5 * ones(Nsim+1, 1);

fprintf('Running Linear MPC...\n');
[x_lmpc, u_lmpc] = runLinearMPC(params, Ts, Nsim, theta_l_ref, umin, umax, dumax);

fprintf('Running Nonlinear MPC...\n');
[x_nmpc, u_nmpc] = runNMPC(params, Ts, Nsim, theta_l_ref, umin, umax, dumax);

metrics_lmpc = computeMetrics(x_lmpc, u_lmpc, theta_l_ref, params, Ts);
metrics_nmpc = computeMetrics(x_nmpc, u_nmpc, theta_l_ref, params, Ts);

fprintf('                         PERFORMANCE COMPARISON\n');
fprintf('%-25s | %-15s | %-15s\n', 'Metric', 'Linear MPC', 'NMPC');
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%-25s | %-15.4f | %-15.4f\n', 'RMS Load Error (rad)', metrics_lmpc.rms_error, metrics_nmpc.rms_error);
fprintf('%-25s | %-15.4f | %-15.4f\n', 'Steady-State Error (rad)', metrics_lmpc.ss_error, metrics_nmpc.ss_error);
fprintf('%-25s | %-15.3f | %-15.3f\n', 'Settling Time (s)', metrics_lmpc.settling_time, metrics_nmpc.settling_time);
fprintf('%-25s | %-15.1f | %-15.1f\n', 'Time Engaged (%%)', metrics_lmpc.time_engaged*100, metrics_nmpc.time_engaged*100);
fprintf('%-25s | %-15.2f | %-15.2f\n', 'Overshoot (%%)', metrics_lmpc.overshoot, metrics_nmpc.overshoot);
fprintf('%-25s | %-15.4f | %-15.4f\n', 'Control Effort', metrics_lmpc.control_effort, metrics_nmpc.control_effort);
fprintf('%-25s | %-15.4f | %-15.4f\n', 'Final Position (rad)', x_lmpc(3,end), x_nmpc(3,end));


fprintf('\n*** WINNER FOR EACH METRIC ***\n');
if metrics_lmpc.rms_error < metrics_nmpc.rms_error
    fprintf('  RMS Error:        Linear MPC\n');
else
    fprintf('  RMS Error:        NMPC\n');
end
if abs(metrics_lmpc.ss_error) < abs(metrics_nmpc.ss_error)
    fprintf('  Steady-State:     Linear MPC\n');
else
    fprintf('  Steady-State:     NMPC\n');
end
if metrics_lmpc.settling_time < metrics_nmpc.settling_time
    fprintf('  Settling Time:    Linear MPC\n');
else
    fprintf('  Settling Time:    NMPC\n');
end
if metrics_lmpc.IAE < metrics_nmpc.IAE
    fprintf('  IAE:              Linear MPC\n');
else
    fprintf('  IAE:              NMPC\n');
end


figure('Name', 'Linear MPC vs NMPC Comparison', 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);

c_lmpc = [0.8 0.2 0.2];
c_nmpc = [0.2 0.2 0.8];

subplot(3,2,1);
plot(t, theta_l_ref, 'k--', 'LineWidth', 2); hold on;
plot(t, x_lmpc(3,:), 'Color', c_lmpc, 'LineWidth', 1.5);
plot(t, x_nmpc(3,:), 'Color', c_nmpc, 'LineWidth', 1.5);
grid on; ylabel('\theta_l (rad)');
legend('Reference', 'Linear MPC', 'NMPC', 'Location', 'southeast');
title('Load Position Tracking');

subplot(3,2,2);
plot(t, x_lmpc(1,:), 'Color', c_lmpc, 'LineWidth', 1.5); hold on;
plot(t, x_nmpc(1,:), 'Color', c_nmpc, 'LineWidth', 1.5);
grid on; ylabel('\theta_m (rad)');
legend('Linear MPC', 'NMPC');
title('Motor Position');

subplot(3,2,3);
stairs(t(1:end-1), u_lmpc, 'Color', c_lmpc, 'LineWidth', 1.2); hold on;
stairs(t(1:end-1), u_nmpc, 'Color', c_nmpc, 'LineWidth', 1.2);
yline(umin, 'k--'); yline(umax, 'k--');
grid on; ylabel('u (N·m)');
legend('Linear MPC', 'NMPC');
title('Control Input');

subplot(3,2,4);
gap_lmpc = x_lmpc(1,:) - x_lmpc(3,:);
gap_nmpc = x_nmpc(1,:) - x_nmpc(3,:);
plot(t, gap_lmpc, 'Color', c_lmpc, 'LineWidth', 1.2); hold on;
plot(t, gap_nmpc, 'Color', c_nmpc, 'LineWidth', 1.2);
yline(params.d, 'r--', 'LineWidth', 1.5);
yline(-params.d, 'r--', 'LineWidth', 1.5);
patch([t(1) t(end) t(end) t(1)], [-params.d -params.d params.d params.d], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on; ylabel('Gap (rad)');
legend('Linear MPC', 'NMPC');
title('Relative Displacement (Gap)');

subplot(3,2,5);
err_lmpc = x_lmpc(3,:)' - theta_l_ref;
err_nmpc = x_nmpc(3,:)' - theta_l_ref;
plot(t, err_lmpc, 'Color', c_lmpc, 'LineWidth', 1.2); hold on;
plot(t, err_nmpc, 'Color', c_nmpc, 'LineWidth', 1.2);
yline(0, 'k--');
grid on; ylabel('Error (rad)'); xlabel('Time (s)');
legend('Linear MPC', 'NMPC');
title('Tracking Error');


sgtitle('Comparison: Linear MPC vs Nonlinear MPC', 'FontSize', 14, 'FontWeight', 'bold');

function [x_hist, u_hist] = runLinearMPC(params, Ts, Nsim, theta_l_ref, umin, umax, dumax)
    Jm = params.Jm; Jl = params.Jl;
    Bm = params.Bm; Bl = params.Bl;
    Ks = params.Ks; Bs = params.Bs;
    d = params.d;
    
    A = [0, 1, 0, 0;
        -Ks/Jm, -(Bm+Bs)/Jm, Ks/Jm, Bs/Jm;
         0, 0, 0, 1;
         Ks/Jl, Bs/Jl, -Ks/Jl, -(Bl+Bs)/Jl];
    B = [0; 1/Jm; 0; 0];
    C = [eye(4); 1, 0, -1, 0];
    D = zeros(5, 1);
    
    sysc = ss(A, B, C, D);
    sysd = c2d(sysc, Ts, 'zoh');
    
    mpcobj = mpc(sysd, Ts, 50, 10);
    
    mpcobj.Weights.OV = [50 1 100 1 300];
    mpcobj.Weights.MV = 0;
    mpcobj.Weights.MVRate = 0.05;
    
    mpcobj.OV(5).Min = -10;
    mpcobj.OV(5).Max = 10;
    mpcobj.MV.Min = umin;
    mpcobj.MV.Max = umax;
    mpcobj.MV.RateMin = -dumax;
    mpcobj.MV.RateMax = dumax;
    
    gap_ref = d + 0.005;
    yref = zeros(Nsim+1, 5);
    for k = 1:Nsim+1
        yref(k,:) = [0, 0, theta_l_ref(k), 0, gap_ref];
    end
    
    x_hist = zeros(4, Nsim+1);
    u_hist = zeros(Nsim, 1);
    x_hist(:,1) = [0;0;0;0];
    
    xmpc = mpcstate(mpcobj);
    
    n_sub = 10;
    dt = Ts/n_sub;
    
    for k = 1:Nsim
        y_meas = C * x_hist(:,k);
        u = mpcmove(mpcobj, xmpc, y_meas, yref(k,:)');
        u = max(min(u, umax), umin);
        u_hist(k) = u;
        
        x_temp = x_hist(:,k);
        for s = 1:n_sub
            x_temp = nonlinearDynamics(x_temp, u, dt, params);
        end
        x_hist(:,k+1) = x_temp;
    end
end

function [x_hist, u_hist] = runNMPC(params, Ts, Nsim, theta_l_ref, umin, umax, dumax)
    Np = 30;
    Nc = 10;
    
    Q = diag([10, 0.1, 200, 1]);
    R = 0.001;
    Rd = 0.01;
    Q_dz = 50;
    
    x_hist = zeros(4, Nsim+1);
    u_hist = zeros(Nsim, 1);
    x_hist(:,1) = [0;0;0;0];
    
    u_prev = 0;
    u_warm = 0.5 * ones(Nc, 1);
    
    n_sub = 10;
    dt = Ts/n_sub;
    
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
                           'MaxIterations', 100, 'OptimalityTolerance', 1e-4);
    
    for k = 1:Nsim
        x_current = x_hist(:,k);
        
        ref_horizon = zeros(4, Np);
        for j = 1:Np
            idx = min(k+j, Nsim+1);
            ref_horizon(:,j) = [theta_l_ref(idx) + params.d*1.5; 0; theta_l_ref(idx); 0];
        end
        
        costFun = @(u_seq) nmpcCost(u_seq, x_current, u_prev, ref_horizon, params, Ts, Np, Nc, Q, R, Rd, Q_dz);
        
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
            fprintf('  NMPC Step %d/%d\n', k, Nsim);
        end
    end
end

function J = nmpcCost(u_seq, x0, u_prev, ref_horizon, params, Ts, Np, Nc, Q, R, Rd, Q_dz)
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

function metrics = computeMetrics(x_hist, u_hist, theta_l_ref, params, Ts)
    d = params.d;
    Nsim = length(u_hist);
    
    gap = x_hist(1,:) - x_hist(3,:);
    in_gap = abs(gap) <= d;
    
    load_error = x_hist(3,:)' - theta_l_ref;
    metrics.rms_error = sqrt(mean(load_error.^2));
    metrics.ss_error = mean(load_error(end-100:end));
    
    final_ref = theta_l_ref(end);
    settled = abs(x_hist(3,:) - final_ref) < 0.02*abs(final_ref);
    idx = find(settled, 1, 'first');
    if isempty(idx)
        metrics.settling_time = Inf;
    else
        metrics.settling_time = idx * Ts;
    end
    
    metrics.time_engaged = sum(~in_gap) / length(in_gap);
    metrics.overshoot = max(0, (max(x_hist(3,:)) - final_ref)/final_ref*100);
    metrics.control_effort = sum(u_hist.^2) * Ts;
    metrics.IAE = sum(abs(load_error)) * Ts;
end
