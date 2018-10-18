#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}
// Predict the State
void KalmanFilter::Predict() {
	// Combine Knowledge of the World (F) with
  	// Prior State Info (x) in order to get
  	// New State Info
  	x_ = F_*x_;
  	// Update the State Uncertainty (P) by 
    // Using Knowledge of the World (F), the
    // Process Noise, and the Prior State Uncertainty
  	P_ = F_*P_*F_.transpose() + Q_;
}
// Update the state by using Kalman Filter Equations
void KalmanFilter::Update(const VectorXd &z) {  
  	// KF Equations
  	VectorXd zed = H_*x_;
	VectorXd y = z - zed;
  	UpdateWithY(y);
}

void KalmanFilter::UpdateWithY(const VectorXd &y) {
 	MatrixXd S = H_*P_*H_.transpose() + R_;
  	MatrixXd K =  P_*H_.transpose()*S.inverse();
  	x_ = x_ + K*y;
  
  	// Init Identity Matrix Based on Shape of x_
    long n = x_.size();
  	MatrixXd I = MatrixXd::Identity(n, n);
  	P_ = (I - K*H_)*P_; 
}

// Update the state by using Extended Kalman Filter equations
void KalmanFilter::UpdateEKF(const VectorXd &z) {
	// Precalculations
  	float px = x_(0), 
          py = x_(1), 
          vx = x_(2), 
          vy = x_(3);
  	// Convert from Cartesian to Polar Coordinates
  	float r = sqrt(px*px + py*py);
  	float th = atan2(py, px);
  	float r_dot = 0;
  
  	// Check if rho is 0 (or close to it)
  	if(fabs(r) >= 0.0001) {
      r_dot = (px*vx + py*vy)/r;
    }
  	// Package conversions into a vector and Update KF
  	VectorXd zed(3);
  	zed << r, th, r_dot;
  	VectorXd y = z - zed;
  	// Update KF
  	UpdateWithY(y);
}
