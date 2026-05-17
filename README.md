# MAE 506 - LMPC vs NMPC Robot Arm Control

This repository contains MATLAB implementations and project artifacts for comparing LMPC and NMPC controllers for a robot arm.

## Expected Project Files
- `LMPC_vs_NMPC_Final.m`
- `MAE_LTIMPC.m`
- `NMPC_DeadzonePenalty.m`
- `NMPC_Optimized.m`
- `NMPC_RobotArm_Animation (1).mp4`
- `MAE 506 Project report Team 12.docx.pdf`

# MAE 506 - LMPC vs NMPC Robot Arm Control

This project compares Linear MPC (LMPC) and Nonlinear MPC (NMPC) for a two-mass robot arm model with backlash (dead-zone) in the shaft coupling.  
The goal is to track a load position reference while handling actuator limits and reducing tracking error.

## Problem Setup
The model represents:
- Motor side inertia and damping
- Load side inertia and damping
- Elastic shaft with damping
- Dead-zone/backlash behavior using gap threshold `d`

State vector used in the scripts:
`x = [theta_m; omega_m; theta_l; omega_l]`

Where:
- `theta_m`, `omega_m`: motor position and velocity
- `theta_l`, `omega_l`: load position and velocity

## What Each Script Does

### `MAE_LTIMPC.m`
Implements a linear MPC controller using a linearized engaged model (`ss`, `c2d`, `mpc`).  
It simulates the *nonlinear* plant in closed loop and reports:
- RMS tracking error
- Time spent in/out of dead-zone
- Gap behavior and control effort

Use this as a baseline LMPC result.

### `NMPC_DeadzonePenalty.m`
Implements NMPC with `fmincon` and a dedicated dead-zone penalty term in the cost function.  
The penalty increases when the motor-load gap is inside the backlash region, pushing behavior toward engagement.

Good for testing explicit dead-zone shaping in the objective.

### `NMPC_Optimized.m`
An improved NMPC configuration with:
- Longer prediction horizon / tuned horizons
- Input-rate inequality constraints
- Terminal cost weighting
- Warm-started optimization

Use this when you want better final tracking quality and smoother control transitions.

### `LMPC_vs_NMPC_Final.m`
Runs both LMPC and NMPC pipelines and prints side-by-side metrics:
- RMS load error
- Steady-state error
- Settling time
- Overshoot
- Control effort

Also generates comparison plots for positions, input torque, gap, and error.

### `NMPC_RobotArm_Animation (1).mp4`
Visualization video of robot arm behavior under controller action.

### `MAE 506 Project report Team 12.docx.pdf`
Final report describing approach, implementation, and outcomes.

## MATLAB Requirements
- MATLAB (R2021a+ recommended)
- Model Predictive Control Toolbox (`mpc`, `mpcmove`)
- Optimization Toolbox (`fmincon`)
- Control System Toolbox (`ss`, `c2d`)

## Recommended Run Order
1. Run `MAE_LTIMPC.m` for baseline linear MPC behavior.
2. Run `NMPC_DeadzonePenalty.m` to inspect dead-zone-aware NMPC.
3. Run `NMPC_Optimized.m` for tuned NMPC performance.
4. Run `LMPC_vs_NMPC_Final.m` for direct comparison tables and plots.

## Outputs You Should Expect
- Console metrics (error, settling time, overshoot, engagement percentage)
- Time-series plots for:
  - Motor/load positions
  - Control torque
  - Backlash gap (`theta_m - theta_l`)
  - Tracking error

