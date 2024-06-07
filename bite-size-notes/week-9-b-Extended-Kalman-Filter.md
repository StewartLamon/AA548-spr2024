<h1>Extended Kalman Filter</h1>

<h3>Objectives</h3>

The objectives of these notes is to highlight the important differences when applying an extended kalman filter as compared to the standard linear kalman filter.

<h3>Introduction</h3>

The Extended Kalman Filter (EKF) is a powerful extension of the Kalman Filter designed to handle nonlinear systems. While the standard Kalman Filter excels in linear scenarios, many real-world systems exhibit nonlinear behavior that the standard filter cannot address adequately. The EKF overcomes this limitation by linearizing the nonlinear system around the current estimate using a first-order Taylor series expansion at each time step.

<h3>Preliminaries: Notation and problem setup</h3>

We will assume Nonlinear Discrete Dynamics with zero mean gaussian process and measurement noise

$$\begin{align*}
x_{k+1}&=f(x_k,u_k,w_k) &&\qquad w_k\sim N(0,Q)\\
y_{k+1}&=g(x_k,v_k) &&\qquad v_k\sim N(0,R)
\end{align*}$$

where f() describes how the dynamics propagate forward in time and g() describes how the states map to measurements of the system.

At each step we will linearize and use the standard kalman filter equations. To linearize we will take the jacobian of both f() and g() w.r.t. the state x, and noise w/v.

$$\begin{align*}
F_{k}^{x}=\frac{\partial{f}}{\partial{x}}\bigg{|}\_{x_k,u_k,w_k}&
F_{k}^{w}=\frac{\partial{f}}{\partial{w}}\bigg{|}\_{x_k,u_k,w_k}&
G_{k}^{x}=\frac{\partial{g}}{\partial{x}}\bigg{|}\_{x_k,u_k,v_k}&
G_{k}^{v}=\frac{\partial{g}}{\partial{v}}\bigg{|}\_{x_k,u_k,v_k}
\end{align*}$$

Here subscript denotes the timestep and the superscript denotes what the jacobian is being taken with respect to. Additionally if T is in the superscript like $G_{k}^{xT}$, then that means it is the transpose of $G_{k}^{x}$

For generalization it is assumed that the noise is not additive. If the noise is additive then the partials with respect to w or v will be identity matrices.

<h3>General Extended Kalman Filter Algorithm</h3>

Predict:

$$\begin{align*}
\mu_{k+1}^p&=f(\mu_k,u_k)\\
\Sigma_{k+1}^p&=F_k^x\Sigma_k F_k^{xT}+F_k^wQF_k^{wT}\\
y_{k+1}^p&=g(\mu_{k+1}^p)
\end{align*}$$

Update:

$$\begin{align*}
K_{k+1}&=\Sigma_{k+1}^pG_{k+1}^{xT}(G_{k+1}^{x}\Sigma_{k+1}^p G_{k+1}^{xT}+G_{k+1}^{v}RG_{k+1}^{vT})^{-1}\\
\mu_{k+1}&=\mu_{k+1}^p+K_{k+1}(y_{k+1}-y_{k+1}^p)\\
\Sigma_{k+1}&=\Sigma_{k+1}^p-K_{k+1}(G_{k+1}^x\Sigma_{k+1}^pG_{k+1}^{xT}+G_{k+1}^vRG_{k+1}^{vT})K_{k+1}^T
\end{align*}$$

Notice that we use the full nonlinear system in the mean prediction step but we have to use the linearized system when looking at how the covariance propagates through the dynamics. This can lead to errors if the dynamics are not well characterized for small pertubations by a linear system around that point. So it's important to make sure your system can be linearized well around expected states.

Simple example: Problem Setup for Nonlinear Pendulum

The nonlinear dynamics of the pendulum can be written as

$$\frac{d}{dt}\begin{bmatrix}
x\\
\dot{x}\end{bmatrix}
=\begin{bmatrix}
\dot{x} \\
-\sin(x)
\end{bmatrix}
+
\begin{bmatrix}
0 \\
u
\end{bmatrix}
+
\begin{bmatrix}
0 \\
w
\end{bmatrix}$$

Which can be discretized using a zero-order hold to

$$\begin{align*}
\begin{bmatrix}
x_{k+1} \\
\dot{x}_{k+1}
\end{bmatrix}
&=\begin{bmatrix}  
x_k+\dot{x}_k\Delta t \\
\dot{x}_k-\sin(x_k)\Delta t
\end{bmatrix}
+\begin{bmatrix}  
0 \\
u\Delta t
\end{bmatrix}
+\begin{bmatrix}
0 \\
w\Delta t
\end{bmatrix} \\
y\_{k+1}&=x\_{k+1}+v\_{k+1}
\end{align*}$$

Notice that we are only measuring the angle of the pendulum (x)

we can calculate our jacobians as

$$\begin{align*}
F_k^{x}&=\frac{\partial f}{\partial x}=\begin{bmatrix}1&\Delta{t}\\ 
-cos(x_k)\Delta{t}&{1}
\end{bmatrix}\\
F^{w}&=\frac{\partial f}{\partial w} = I\\
G^{x}&=\frac{\partial g}{\partial x} = \begin{bmatrix} 1 \\ 0 \end{bmatrix}\\
G^{v}&=\frac{\partial g}{\partial v} = I
\end{align*}$$

In this case we can see that our dynamics f(x) are nonlinear but our measurement function g(x) is linear. What has resulted is that $G^{x}$ is constant whereas $F_k^x$ will change depending on the state $x_k$, hence it has the subscript k. Also, the noise was assumed to be additive and therefore the jacobians w.r.t. the noise were Identity matrices.

This was simulated at two different noise levels on both a kalman filter linearized around the fixed point and an extended kalman filter\\
First with a low noise level $R=0.015$

![image](https://github.com/StewartLamon/AA548-spr2024/assets/128524152/14fed07f-cdee-4b08-9489-6ee7ac2babb7)
![image](https://github.com/StewartLamon/AA548-spr2024/assets/128524152/01430859-79ff-404d-a1eb-f1d2db4853e3)

The linear kalman filter performs pretty well. Let's look at it again with increased noise $R=0.15$.

![image](https://github.com/StewartLamon/AA548-spr2024/assets/128524152/5e530afe-90aa-4848-96d4-d4e298956a58)
![image](https://github.com/StewartLamon/AA548-spr2024/assets/128524152/1d51f0e6-f938-4ebd-86ea-5a39b2971e96)

The linear kalman filter becomes unstable! It cannot accurately estimate how the covariance matrix passes through the dynamics.

<h3>On optimality and stability:</h3>

One important consequence of linearizing is that this is no longer a truely optimal filter since we are no longer using the actual dynamics but an approximation during the covariance prediction/update steps. This is very important to take into account because the extended kalman filter cannot be used if the nonlinear dynamics is highly nonlinear, such that a local linearization would be invalid and give erroneous estimates.

<h3>Conclusion</h3>

The Extended Kalman Filter (EKF) is a powerful extension of the Linear Kalman Filter, capable of handling non-linear dynamics and measurements. By approximating the change in covariance through linearization of the non-linear functions f() and g() at each time step using the current estimated state, the EKF effectively manages the complexities introduced by non-linearity. This approach ensures that the filter remains computationally feasible while providing accurate state estimates.

<h3>Refences</h3>

[1] Extended Kalman Filters
- MATLAB & Simulink. (n.d.). https://www.mathworks.com/help/fusion/ug/extended-kalman-filters.html


