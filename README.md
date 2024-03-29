# Attitude estimation

This package provides attitude and rate estimation algorithms for aerial robotics applications. 
Observations are gyro-free and are based on two vectorial measurements.
Both the Extended Kalman filter and the Unscented Kalman filter are modified to include a symplectic integration step. This is how momentum and energy are preserved naturally.
A geometric filter and a second-order-optimal minimum energy filter are compared against the aforementioned Gaussian-approximate filters.

For both case studies, the measurement noise and model uncertainty are initially expressed under the stochastic reasoning. Subsequently, to stress the significance of the dual optimal control formulation, we replace the model error with an unknown deterministic disturbance that exerts the actual system.

# Motivation
Stochastic modeling appears to be problematic for this specific problem. In contrast to the Kalman-based filters, deterministic attitude and rate filtering can be easily formulated in a coordinate-free fashion as a dual optimal control problem. This approach provides the ability to overcome the well-known singularity issues imposed by the interaction of the pre-fabricated Bayesian architectures with the state constraints.

# Results
Both deterministic filters outperform the Gaussian-approximate filters. This is mainly due to the fact that the necessary re-projection step induces a bias in the orientation estimates. Regarding the stochastic filters we conclude that due to the computational overhead of the UKF, the simplicity of the Jacobian matrix calculations, and the quasi-linear nature of the quaternion dynamics, the EKF is 
a better choice for the task.
 
Attitude and rate estimation errors for UAV: Stochastic modeling (left), Deterministic modeling (right).

<p float="left">
  <img src="figures_png/orientation.png" width="400" height="220"/>
  <img src="figures_png/orientationmodel.png" width="400" height="220"/> 
</p>

<p float="left">
  <img src="figures_png/WX.png" width="400" height="220"/>
  <img src="figures_png/WXmodel.png" width="400" height="220"/> 
</p>

<p float="left">
  <img src="figures_png/WY.png" width="400" height="220" style="filter: brightness(0.1);"/>
  <img src="figures_png/WYmodel.png" width="400" height="220"/> 
</p>

<p float="left">
  <img src="figures_png/WZ.png" width="400" height="220"/>
  <img src="figures_png/WZmodel.png" width="400" height="220"/> 
</p>

