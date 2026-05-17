
clear; close all; clc;

%% Physical parameters
Jm = 0.02;    % motor inertia [kg*m^2]
Jl = 0.05;    % load inertia [kg*m^2]
Bm = 0.01;   % motor viscous damping
Bl = 0.01;   % load viscous damping
Ks = 5;      % shaft stiffness (when engaged)
Bs = 0.05;    % shaft damping (when engaged)
d  = 0.03;    % dead-zone half-width (rad)  

%% Linearized continuous-time state-space (engaged)
% states x = [theta_m; theta_m_dot; theta_l; theta_l_dot]
A = [ 0,              1,           0,              0;
     -Ks/Jm,  -(Bm+Bs)/Jm,    Ks/Jm,        Bs/Jm;
      0,              0,           0,              1;
      Ks/Jl,      Bs/Jl,     -Ks/Jl,     -(Bl+Bs)/Jl ];
B = [0; 1/Jm; 0; 0];
C_states = eye(4);
C_gap = [1, 0, -1, 0];        % single row
C = [C_states; C_gap];        % 5 outputs (4 states + gap)
D = zeros(size(C,1), size(B,2));

%% Discretize plant 
Ts = 0.01;    % sample time(s)
sysc = ss(A,B,C,D);
sysd = c2d(sysc, Ts, 'zoh');
Ad = sysd.A; Bd = sysd.B; Cd = sysd.C; Dd = sysd.D;

%% Create MPC controller
PredictionHorizon = 50;   % Np
ControlHorizon    = 10;    % Nc

mpcobj = mpc(sysd, Ts, PredictionHorizon, ControlHorizon);

% Weights:

Q_states = [50 1 100 1];   % weight for [theta_m, theta_m_dot, theta_l, theta_l_dot]
w_gap   = 300;              % weight for the gap output
mpcobj.Weights.OV = [Q_states, w_gap];

% Penalize MV and MVRate: R on Delta u is implemented via MVRate weight
mpcobj.Weights.MV = 0;         % small direct MV penalty
mpcobj.Weights.MVRate = 0.05;   % penalize changes in u 

% Soft constraint ECR weight 
mpcobj.OV(5).Min = -10; % we'll use references rather than hard OV min/max
mpcobj.OV(5).Max = 10;

% Actuator constraints
umin = -1;   % torque min 
umax =  1;   % torque max 
dumax = 10;  % max rate change

mpcobj.MV.Min = umin;
mpcobj.MV.Max = umax;

mpcobj.MV.RateMin = -dumax*Ts;
mpcobj.MV.RateMax =  dumax*Ts;


%% Simulation parameters
Tsim = 6;              % seconds
Nsim = round(Tsim/Ts);
t = (0:Nsim)'*Ts;

% Reference trajectories: aim to track a load position theta_l_ref
theta_l_ref_traj = 0.5*ones(Nsim+1,1);   % step reference for load (rad)
gap_margin = 0.005;  % margin beyond dead-zone
gap_ref_value = sign(1)*(d + gap_margin); % target gap (pick +d or -d depending on initial)
% Build reference for entire output vector (OV = [theta_m; theta_m_dot; theta_l; theta_l_dot; gap])
yref = zeros(Nsim+1, size(C,1));
for k=1:Nsim+1
    yref(k,1) = 0;                     % theta_m ref 
    yref(k,2) = 0;                     % theta_m_dot ref
    yref(k,3) = theta_l_ref_traj(k);   % theta_l ref 
    yref(k,4) = 0;                     % theta_l_dot ref
    yref(k,5) = gap_ref_value;         % the gap to be >= d 
end

%% Closed-loop simulation on the FULL nonlinear model 
% Nonlinear piecewise model function (engaged vs disengaged)
plant_step = @(x,u) deal( Ad*x + Bd*u ); % we'll simulate discrete-time nonlin model manually below

% Jm*theta_m_dd = u - Bm*theta_m_d - tau_s
% Jl*theta_l_dd = tau_s - Bl*theta_l_d
% tau_s = 0 if |theta_m-theta_l| <= d
% tau_s = Ks*(theta_m-theta_l) + Bs*(theta_m_d - theta_l_d) if outside dead-zone

% We'll integrate with simple Euler at Ts (fine if Ts small); for better accuracy use ode45.

% Initialize states (continuous-state vector x_c)
x = zeros(4, Nsim+1);
% initial conditions 
x(:,1) = [0; 0; 0; 0];   

% mpc state object for internal controller state
xmpc = mpcstate(mpcobj);

% history store
u_hist = zeros(Nsim,1);
y_hist = zeros(Nsim+1, size(C,1));
y_hist(1,:) = (C * x(:,1))';

for k = 1:Nsim
    % form current y measurement 
    y_meas = (C * x(:,k));
    
    % current reference for horizon
    refk = yref(k,:)';   % mpcmove accepts current-time reference 
    
    % compute control action with mpcmove
    mv = mpcmove(mpcobj, xmpc, y_meas, refk);
    u = mv;  % scalar torque command
    
    % apply MV saturation to be safe
    u = max(min(u, umax), umin);
    
    % integrate continuous nonlinear dynamics over one sample time Ts using simple Euler:
    th_m = x(1,k); thm_d = x(2,k);
    th_l = x(3,k); thl_d = x(4,k);
    rel  = th_m - th_l;
    reld = thm_d - thl_d;
    
    if abs(rel) <= d
        tau_s = 0;
    else
        tau_s = Ks*rel + Bs*reld;
    end
    
    % continuous derivatives
    thm_dd = (u - Bm*thm_d - tau_s) / Jm;
    thl_dd = (tau_s - Bl*thl_d) / Jl;
    
    % Euler integration (simple)
    x(1,k+1) = th_m + Ts*thm_d;
    x(2,k+1) = thm_d + Ts*thm_dd;
    x(3,k+1) = th_l + Ts*thl_d;
    x(4,k+1) = thl_d + Ts*thl_dd;
    
    % store
    u_hist(k) = u;
    y_hist(k+1,:) = (C * x(:,k+1))';
end

%% Metrics 
gap = y_hist(:,5);            % this is theta_m - theta_l
in_gap = abs(gap) <= d;
time_in_gap_fraction = sum(in_gap)/length(in_gap);
rms_gap = sqrt(mean(gap.^2));
load_tracking_error = x(3,:)' - theta_l_ref_traj; % theta_l - ref
rms_load_error = sqrt(mean(load_tracking_error.^2));

fprintf('Time fraction in gap (|g|<=d): %.3f\n', time_in_gap_fraction);
fprintf('RMS gap: %.4f rad\n', rms_gap);
fprintf('RMS load tracking error: %.4f rad\n', rms_load_error);


%% Plots 
figure('Name','MPC Closed-loop Simulation','Units','normalized','Position',[0.1 0.1 0.6 0.75]);

subplot(4,1,1);
plot(t, x(1,:), 'LineWidth', 1.5); hold on;
plot(t, x(3,:), '--', 'LineWidth', 1.5);
grid on;
ylabel('\theta (rad)');
legend('\theta_m','\theta_l','Location','best');
title('Positions');

subplot(4,1,2);
stairs(t(1:end-1), u_hist, 'LineWidth', 1.5); grid on;
ylabel('u (N\cdotm)');
title('Control input (torque)');

subplot(4,1,3);
plot(t, gap, 'LineWidth', 1.5); hold on;
h1 = yline(d, 'r--');
h2 = yline(-d, 'r--');
ax = gca;
x_text = ax.XLim(2) - 0.05*(ax.XLim(2)-ax.XLim(1));
text(x_text, d, '  +d', 'Color','r','FontWeight','bold', 'HorizontalAlignment','left');
text(x_text, -d, '  -d', 'Color','r','FontWeight','bold', 'HorizontalAlignment','left');
grid on;
ylabel('gap = \theta_m - \theta_l (rad)');
title('Relative displacement (gap)');

subplot(4,1,4);
plot(t, double(in_gap), 'LineWidth', 1.5); grid on;
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'out','in'});
ylabel('gap');
xlabel('time (s)');

sgtitle('MPC Closed-loop Simulation (Linear MPC on engaged linearization)');

%% End of script
