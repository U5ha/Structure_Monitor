# Structure_Monitor
Monitoring Structural health based on Vibration data of buildings and other structures.

**System States Tracking**
The Kalman Filter is an efficient optimal estimator (a set of mathematical equations) that provides a recursive computational methodology for estimating the state of a discrete-data controlled process from measurements that are typically noisy, while providing an estimate of the
uncertainty of the estimates (Thomson and Emery (2014)). 
Kalman Filter is widely used nowadays in the structural health monitoring field, where the system response for the kth step can be estimated using the data upto (k − 1)th step. This estimated data is updated after the kth measurement by giving proper weightages to the estimated and measured data.
• The initial guess of the state vector and the covariance is taken according to the question.
• The estimation of the state estimate and covariance matrix, for the next step is done using the data upto previous step.
• Measurement of the data is done for the current step.
• Kalman gain is calculated.
• The state estimate and the covariance matrix is updated by giving proper weightages to the previously estimated and the measured values.
It is seen that in case of velocity measurements the predicted and measured value have similar shape but the predicted and measured value are shifted apart. The reason behind this is the assumption of the initial state vector. Farther the initial assumption of the velocity value, more the accumulated error, because displacement is integration of velocity, so the error keeps on increasing. On the other hand, when the calculations are done by measuring the displacement
data, the error keeps on decreasing from displacement to velocity to acceleration as the latter is obtained by differentiating the former.

**Process Noise Effects**
The system can be represented by the continuous state-space model by using Z˙ and Z, which can be simplified in the form and the matrix F is formulated. All the formulations given above are for the continuous state-space model. To convert this model to discrete state space model, either the inbuilt (MATLAB (2022)) function can be used or the formulation given could be used.
Interpretation of the measurement data: Two cases of measurements are given in the problem statement, namely, acceleration and displacement measurements.
 The variance for displacement vector can easily be obtained from the measured data, but, for the variance in velocity data, the average value of variances for the displacement and acceleration is taken. The covariance matrix will thus consist of the various variances. Initially, the process noise Q is assumed to be zero, and the model is assumed to be perfect.
The applied ground motion and the measurement data are of different frequencies, so, the ground motion data is interpolated to obtain equal number of datasets for both ground motion and measured data. The states are tracked using Kalman filter from the measured data for accelerations and displacements. The process noise depends on the error in modelling of the state-space system. The noise matrix is a diagonal matrix with noises in all the states considered to be equal.

**Extended Kalman Filter**
The use of Kalman filter can be extended to estimation of the state from the noisy measurements when the structural parameters are unknown. The only change that occurs in the state vector in this case is the inclusion of structural parameters, and the new state vector is called augmented state vector. This means that the structural parameters will also be estimated at each step of calculation. This concept of augmented state vector is used in this assignment to track the states
of a dynamic system subjected to external force.
The unknown parameters are to be estimated using the noisy measurements and system is to be tracked using Extended Kalman filter. The frequency of the measurement data is 200 Hz and that of the applied ground excitation is 50 Hz, so the applied ground excitation is to be interpolated for the intermediate time steps.
The point to be noted here is that due to the unknown structural parameters the system becomes non-linear . This is the continuous system formulation which is to be converted to discrete system in order to compute this in MATLAB.The discretization is done by second order Runge-Kutta method using Ralston’s method, the formulation of F from f using Ralston’s method is done in MATLAB using symbolic substitution.

**Unscented Kalman Filter**
The formulation of state space equations is same as Extended Kalman Filter. This means that the functions F and H remains unchanged. Now, for the calculation of the functional values of the non-linear function Unscented Kalman Filter is used in this case. Instead of approximating the function till first order using Taylor’s series approximation, the functional values are propagated through the actual non-linear functions by selection 2N + 1 sigma points
where, N is the number of state variables. 
A comparison of the state estimates using EKF and UKF with the original noisy measurements is done. The matlab code contains all the provided choices of P and Q matrices. Each combination of P and Q yields different results for state parameters. When no process noise is assumed the results of stiffness and damping do not converge to the original value. For the same case, the velocities and displacements give results that are haywire. For the estimation problem when
the original values are not known, incorrect assumption of Q will lead to erroneous results.
In case of EKF when Q is chosen between 2 and 3 and state estimates are compared, the displacement and velocity states are close to each other, but, in estimation of parameters we can see that there is large noise for higher value of Q, although the final values are very close.
When very high value of P is chosen, i.e., choice 1, the results of displacements and velocities from UKF are monotonically increasing to very high order which means that the system is unstable. When EKF and UKF are compared the higher values are obtained at the peaks with UKF, but, values are quite close for rest of the values. The convergence of UKF to the original values is faster in comparison to EKF.
