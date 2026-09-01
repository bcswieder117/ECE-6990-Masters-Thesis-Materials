# Disturbance-Aware Model-Based and Model-Free Zero-Sum Pursuit-Evasion Control

## Theory, Algorithm Development, and MATLAB Validation

> **Project Status:** Completed Master of Science project report in Electrical and Computer Engineering at Tennessee Technological University, August 2026. This repository contains the reproducible MATLAB implementations, numerical studies, literature documentation, and supporting materials developed for the report. The report evaluates four disturbance-aware zero-sum control methods in a common MATLAB framework; nonlinear vehicle, QLabs, and physical QCar validation remain future work.

## Overview

This repository documents a theoretical and computational research project involving zero-sum dynamic games, model-based control, model-free control, persistent command disturbances, and pursuit-evasion systems.

The project report investigates how established model-based and model-free zero-sum control algorithms can be extended to account for persistent differences between commanded and realized control inputs.

The work does not replace the original algorithms or claim that policy iteration, value iteration, or Q-learning were developed in this project. Instead, the original solution structures are retained and expanded with:

- persistent command-channel offsets;
- command-offset estimation;
- affine input compensation;
- matched nominal and disturbed simulations;
- closed-loop numerical checks;
- MATLAB figure and data generation.

The model-based methods use known system matrices and generalized algebraic Riccati equation methods. The model-free methods recover the corresponding zero-sum policies from sampled state, action, cost, and successor-state data without directly using the analytical state-transition matrices in the learning update.

The completed report provides a theoretical and MATLAB-based foundation for later nonlinear autonomous-ground-vehicle and QCar studies. Completed physical QCar validation is not claimed.

## Report Structure

The project report is organized around four established zero-sum control methods:

1. Model-based policy iteration.
2. Model-based value iteration.
3. Model-free Q-learning policy iteration.
4. Model-free Q-learning value iteration.

Each method is extended using the same disturbance model, compensation structure, and evaluation criteria.

The objective was to produce matched disturbance-aware implementations of the four algorithms while preserving the mathematical role of each original method.

## Zero-Sum Game Formulation

The nominal discrete-time system is represented as:

```text
x(k+1) = A*x(k) + Bp*uP(k) + Be*uE(k)
```

where:

- `x(k)` is the game state;
- `uP(k)` is the pursuer input;
- `uE(k)` is the evader input;
- `Bp` is the pursuer input matrix;
- `Be` is the evader input matrix.

The pursuer is the minimizing player, while the evader is the maximizing player.

A representative quadratic zero-sum cost is:

```text
J = sum over k of:

    x(k)'*Q*x(k)
  + uP(k)'*RP*uP(k)
  - uE(k)'*RE*uE(k)
```

The feedback policies are written as:

```text
uP(k) = L*x(k)
uE(k) = K*x(k)
```

where:

- `L` is the pursuer feedback gain;
- `K` is the evader feedback gain.

The nominal closed-loop system is:

```text
x(k+1) = (A + Bp*L + Be*K)*x(k)
```

## Persistent Command-Channel Offsets

The disturbance-aware extension assumes that the realized inputs may differ from the commanded inputs.

For the pursuer:

```text
uP_actual(k) = uP_command(k) + bP
```

For the evader:

```text
uE_actual(k) = uE_command(k) + bE
```

The terms `bP` and `bE` represent constant or slowly varying command-channel offsets.

Possible physical interpretations include:

- acceleration-command calibration error;
- steering-center offset;
- persistent actuator mismatch;
- repeatable differences between commanded and realized inputs;
- approximately constant implementation error during an experiment.

Without compensation, the disturbed closed-loop system becomes:

```text
x(k+1) = (A + Bp*L + Be*K)*x(k)
         + Bp*bP
         + Be*bE
```

The feedback gains may remain stabilizing, but the persistent affine term can produce a nonzero steady-state deviation from nominal behavior.

## Command-Offset Estimation

For the model-based methods, the offsets can be estimated from known system matrices and measured state transitions.

The transition residual is:

```text
r(k) = x(k+1)
       - A*x(k)
       - Bp*uP_command(k)
       - Be*uE_command(k)
```

Under the persistent-offset model:

```text
r(k) = Bp*bP + Be*bE + noise(k)
```

The two offsets can be collected into one vector:

```text
bias = [bP; bE]
```

The combined input matrix is:

```text
Bbias = [Bp Be]
```

The calibration relationship is therefore:

```text
r(k) = Bbias*bias + noise(k)
```

A batch of transition residuals can be used to estimate the offsets through least squares.

## Affine Command Compensation

The estimated offsets are denoted by `bP_hat` and `bE_hat`.

The compensated pursuer command is:

```text
uP_command(k) = L*x(k) - bP_hat
```

The compensated evader command is:

```text
uE_command(k) = K*x(k) - bE_hat
```

After the true offsets enter the realized command channels:

```text
uP_actual(k) = L*x(k) - bP_hat + bP
uE_actual(k) = K*x(k) - bE_hat + bE
```

The compensated closed-loop system becomes:

```text
x(k+1) = (A + Bp*L + Be*K)*x(k)
         + Bp*(bP - bP_hat)
         + Be*(bE - bE_hat)
```

When the estimates are exact:

```text
bP_hat = bP
bE_hat = bE
```

the disturbance terms cancel and the nominal closed-loop dynamics are recovered:

```text
x(k+1) = (A + Bp*L + Be*K)*x(k)
```

When the estimates are imperfect, the remaining disturbance depends only on the bias-estimation errors.

## Model-Based Methods

### Model-Based Policy Iteration

The model-based policy-iteration method alternates between policy evaluation and policy improvement.

For fixed pursuer and evader gains, the closed-loop matrix is:

```text
Ai = A + Bp*Li + Be*Ki
```

The policy-evaluation step computes the value matrix associated with the current policies.

The policy-improvement step updates the pursuer and evader gains from the zero-sum saddle-point equations.

The disturbance-aware extension does not replace the Riccati or gain-update equations. It adds:

- persistent command offsets;
- model-based bias calibration;
- affine input compensation;
- nominal and disturbed simulations;
- convergence histories;
- closed-loop stability checks;
- saddle-point curvature checks;
- Riccati residual calculations;
- reproducible MATLAB outputs.

### Model-Based Value Iteration

The model-based value-iteration method repeatedly applies the generalized Riccati value update and then recovers the pursuer and evader gains.

The disturbance-aware value-iteration implementation uses the same:

- command-offset model;
- bias-estimation procedure;
- compensation law;
- simulation cases;
- numerical checks;
- MATLAB output format.

The value-iteration algorithm is treated separately from policy iteration so that the behavior of the two established model-based solution procedures can be compared under the same disturbance conditions.

## Model-Free Methods

The model-free methods recover the zero-sum policies from sampled transition data rather than directly using the analytical system matrices in the Q-learning update.

The joint state-action vector is:

```text
z(k) = [x(k); uP(k); uE(k)]
```

The quadratic Q-function is represented as:

```text
QH(z(k)) = z(k)'*H*z(k)
```

The symmetric matrix `H` contains state, pursuer-input, evader-input, and cross-term information.

A representative block partition is:

```text
H = [Hxx  HxP  HxE
     HPx  HPP  HPE
     HEx  HEP  HEE]
```

The pursuer and evader gains are recovered from the appropriate blocks of `H` using the zero-sum saddle-point conditions.

### Model-Free Q-Learning Policy Iteration

The Q-learning policy-iteration method evaluates the quadratic Q-function associated with fixed pursuer and evader policies using sampled transitions.

It then improves both policies using the estimated Q-function matrix.

The disturbance-aware extension adds:

- commanded input histories;
- realized input histories;
- persistent offset estimation;
- affine input compensation;
- nominal and disturbed simulations;
- feature-rank checks;
- regression-conditioning checks;
- Bellman residuals;
- policy-convergence histories.

### Model-Free Q-Learning Value Iteration

The Q-learning value-iteration method estimates successive quadratic Q-functions from sampled state transitions and updates the pursuer and evader policies from the resulting matrix blocks.

The same disturbance and compensation model is used so that all four algorithms are evaluated under matched assumptions.

## Common Evaluation Cases

Each algorithm is evaluated using the same primary cases.

### 1. Nominal Operation

No command-channel offsets are present:

```text
bP = 0
bE = 0
```

This case establishes the nominal behavior of the algorithm.

### 2. Persistent Offsets Without Compensation

The true offsets enter the realized pursuer and evader inputs, but no estimated bias is subtracted from the commands.

This case measures the effect of persistent command mismatch.

### 3. Persistent Offsets With Estimated Compensation

The estimated offsets are subtracted from the commanded inputs.

This case measures how effectively the estimation and compensation procedure restores nominal behavior.

### 4. Exact Compensation Verification

The true offsets are used in the compensation law.

This is a verification case used to confirm that exact cancellation recovers the nominal trajectory to numerical precision.

## MATLAB Validation

MATLAB is used to examine the theoretical and computational behavior of the extended algorithms.

The completed numerical studies include:

- policy-iteration convergence;
- value-iteration convergence;
- Q-learning regression convergence;
- generalized algebraic Riccati equation residuals;
- closed-loop spectral radius;
- minimizing-player curvature;
- maximizing-player Schur-complement curvature;
- game-matrix conditioning;
- feature-matrix rank;
- regression conditioning;
- Bellman residuals;
- true and estimated command offsets;
- bias-estimation error;
- calibration residuals;
- nominal state trajectories;
- uncompensated state trajectories;
- compensated state trajectories;
- steady-state deviations;
- finite-window zero-sum game cost;
- deviation from nominal behavior;
- exact-compensation verification error.

The MATLAB scripts save reproducible outputs such as:

- MATLAB `.fig` files;
- high-resolution `.png` figures;
- numerical summary `.csv` files;
- complete workspace `.mat` files.

## Completed Four-Method Study

The completed MS project report evaluates all four algorithms under the same plant, game weights, initial state, command offsets, and numerical criteria. The nominal controllers agree to their respective stopping tolerances and share the same closed-loop spectral radius and saddle-point curvature values.

| Method | Outer steps | Defining residual | Closed-loop spectral radius |
| --- | ---: | ---: | ---: |
| Model-based GARE policy iteration | 4 | `6.40e-15` GARE | `0.9819707273` |
| Model-based GARE value iteration | 633 | `9.44e-11` GARE | `0.9819707273` |
| Model-free Q-learning policy iteration | 5 | `4.02e-15` Bellman | `0.9819707273` |
| Model-free Q-learning value iteration | 633 | `3.38e-12` Bellman | `0.9819707273` |

Across all four methods, the uncompensated maximum Q-weighted state deviation is approximately `4.1286e-2`. Known-model transition-residual estimation reduces the maximum deviation to approximately `1.1350e-4` (about 99.725%), while direct command-mismatch estimation reduces it to approximately `4.1668e-6` (about 99.990%). Exact compensation reproduces the nominal trajectory to numerical precision.

The completed study compares:

```text
Nominal operation
Persistent offsets without compensation
Persistent offsets with estimated compensation
Persistent offsets with exact compensation
```

The analytical affine-deviation recurrence agrees with the independently simulated deviations to approximately `1e-15`, confirming that the remaining trajectory error is governed by the residual offset-estimation error rather than by a change in the nominal policy.

## Relationship to Pursuit-Evasion Control

The algorithms are first evaluated on established discrete-time zero-sum state-space systems.

This allows the Riccati, Q-learning, disturbance-estimation, and compensation components to be studied before introducing the additional complexity of nonlinear vehicle dynamics.

The pursuit-evasion interpretation is:

- the pursuer minimizes the game cost;
- the evader maximizes the same cost;
- the pursuer attempts to reduce relative separation;
- the evader attempts to delay or prevent capture;
- persistent command mismatch affects both players;
- compensation is added without replacing the original game solution method.

## Relationship to Autonomous Ground Vehicles

The project report uses autonomous ground-vehicle pursuit-evasion as the motivating application.

A representative vehicle state is:

```text
vehicle_state_i(k) = [
    x_position_i(k)
    y_position_i(k)
    heading_i(k)
    speed_i(k)
]
```

where `i` identifies the pursuer or evader.

A representative vehicle control input is:

```text
vehicle_input_i(k) = [
    acceleration_i(k)
    steering_i(k)
]
```

For a vehicle implementation, the scalar benchmark offsets may be replaced by acceleration and steering offset vectors:

```text
bP = [
    pursuer_acceleration_offset
    pursuer_steering_offset
]
```

```text
bE = [
    evader_acceleration_offset
    evader_steering_offset
]
```

The completed linear zero-sum studies therefore provide the algorithmic and numerical foundation for later nonlinear vehicle extensions.

## Research Focus

- **Zero-sum pursuit-evasion control:** A two-player game involving a minimizing pursuer and a maximizing evader.
- **Model-based policy iteration:** Generalized Riccati policy iteration with command-offset estimation and affine compensation.
- **Model-based value iteration:** Generalized Riccati value iteration with the same disturbance and compensation structure.
- **Model-free Q-learning policy iteration:** Quadratic Q-learning policy iteration using sampled transitions.
- **Model-free Q-learning value iteration:** Quadratic Q-learning value iteration under matched game assumptions.
- **Persistent command disturbances:** Repeatable differences between commanded and realized inputs.
- **Disturbance estimation:** Estimation of command offsets from transition residuals or measured command data.
- **Affine compensation:** Cancellation of estimated command offsets during policy execution.
- **Matched MATLAB validation:** Common simulation cases and numerical criteria across all four algorithms.
- **Vehicle-oriented foundation:** Preparation for later nonlinear kinematic-bicycle, QLabs, or QCar studies.

## Presentation

This work was presented at Tennessee Tech's **2026 Research and Creative Inquiry Symposium**.

[View the Tennessee Tech Research and Creative Inquiry Symposium page](https://www.tntech.edu/research/research-day.php)

The symposium presentation reflects an earlier stage of the MS project. The final August 2026 report places greater emphasis on disturbance-aware extensions of established model-based and model-free zero-sum algorithms and matched MATLAB validation.

## Abstract

This project report investigates disturbance-aware control for a discrete-time, two-player, zero-sum pursuit-evasion game and develops a common framework for extending established model-based and model-free solution methods without replacing their original policy iteration, value iteration, or Q-learning structures. The pursuer is represented as the minimizing player, while the evader is represented as the maximizing player within a linear-quadratic dynamic game. The nominal saddle-point policies are characterized through a generalized algebraic Riccati equation and an equivalent quadratic action-value representation.

Four established methods are considered: model-based generalized Riccati policy and value iteration, as well as model-free quadratic Q-learning policy and value iteration. Each method is extended to include persistent command-channel offsets affecting both the pursuer and evader. Offset estimates are obtained from transition residuals or measured command mismatch, depending on the information structure of the algorithm. Affine compensation is then applied by subtracting the estimated offsets from the commanded control actions.

The resulting residual-bias dynamics directly relate the compensated closed-loop response to the remaining estimation error and establish exact compensation as a verification case. The model-free methods use realized inputs in the sampled transition tuples so that the Bellman regression remains physically consistent with the transition and stage cost produced by those inputs.

MATLAB studies evaluate nominal operation, uncompensated command offsets, estimated-offset compensation, and exact compensation under matched system matrices, game weights, disturbances, initial conditions, and numerical criteria. Validation includes policy and value convergence, generalized Riccati and Bellman residuals, closed-loop spectral radius, saddle-point curvature, matrix conditioning, feature-rank requirements, bias-estimation error, state and input trajectories, and cumulative game cost. The contribution is the unified disturbance-aware extension and matched evaluation of four established zero-sum control methods within a common pursuit-evasion framework. Nonlinear autonomous-ground-vehicle, QLabs, and physical QCar implementations remain future work.

## Scope

The current scope includes:

- one minimizing pursuer and one maximizing evader;
- discrete-time linear or locally linearized systems;
- two-player zero-sum quadratic games;
- model-based policy iteration;
- model-based value iteration;
- model-free Q-learning policy iteration;
- model-free Q-learning value iteration;
- persistent command-channel offsets;
- command-offset estimation;
- affine command compensation;
- MATLAB convergence and closed-loop studies;
- matched numerical evaluation of all four algorithms;
- a foundation for later autonomous-ground-vehicle implementation.

The current scope does not include:

- completed physical QCar experiments;
- completed QLabs deployment;
- a claim that the underlying algorithms were invented in this project;
- deep actor-critic reinforcement-learning methods;
- a full autonomous-driving traffic simulator;
- multi-pursuer or multi-evader coordination;
- global optimality guarantees for nonlinear vehicle dynamics;
- high-fidelity tire, suspension, or three-dimensional vehicle dynamics.

## Research Contribution

The contribution of this project report is not the invention of policy iteration, value iteration, or Q-learning.

The contribution is the systematic extension and evaluation of established model-based and model-free zero-sum algorithms for a disturbance-aware pursuit-evasion problem.

The primary contributions are:

1. A common disturbance-aware zero-sum pursuit-evasion formulation.
2. An extension of model-based policy iteration with persistent command-offset estimation and affine compensation.
3. An extension of model-based value iteration using the same disturbance structure.
4. An extension of model-free Q-learning policy iteration for the same zero-sum game.
5. An extension of model-free Q-learning value iteration under matched disturbance assumptions.
6. A reproducible MATLAB framework for nominal, uncompensated, estimated-compensation, and exact-compensation studies.
7. Matched numerical evidence that all four methods recover the same nominal saddle solution under their respective information structures.
8. A theoretical and computational foundation for later nonlinear vehicle and QCar implementation.

## Future Work

Future work may include:

- nonlinear kinematic-bicycle dynamics;
- dynamic bicycle models;
- vector-valued acceleration and steering offsets;
- online bias estimation;
- slowly time-varying disturbances;
- stochastic command disturbances;
- process and measurement noise;
- partial-state measurements;
- output-feedback formulations;
- actuator saturation;
- communication and command delays;
- multiple pursuers and evaders;
- QLabs simulation;
- physical QCar testing;
- hardware-based command calibration;
- experimental comparison of model-based and model-free implementations.

## Citation

Swieder, Blaine Christopher. *Disturbance-Aware Model-Based and Model-Free Zero-Sum Pursuit-Evasion Control: Theory, Algorithm Development, and MATLAB Validation.* MS project report, Tennessee Technological University, August 2026.
