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
  	// Precalculations
  	auto P_Ht = P_*H_.transpose();
  	auto n = x_.size();
  	auto I = MatrixXd::Identity(n, n);
  
  	// KF Equations
	auto S = H_*P_Ht + R_;
  	auto K = P_Ht*S.inverse();
  	
  	auto K_H = K*H_; // Intermediary Calculation
  	
  	// x_ = x_ + K*(z - H_*x_) -> x_ = x_ + K*z - K*H_*x_ -> x_ = K*z + x_(1 - K*H_)
  	x_ = K*z + x_*(MatrixXd::Constant(1,1,1) - K_H);
  	P_ = (I - K_H)*P_;
}
// Update the state by using Extended Kalman Filter equations
void KalmanFilter::UpdateEKF(const VectorXd &z) {
	// Precalculations
  	auto px = x_(0), 
         py = x_(1), 
         vx = x_(2), 
         vy = x_(3);
  	// Convert from Cartesian to Polar Coordinates
  	auto r = sqrt(px*px + py*py);
  	auto th = atan2(py, px);
  	float r_dot = 0;
  	// Check if rho is 0 (or close to it)
  	if(fabs(r) >= 0.0001) {
      r_dot = (px*vx + py*vy)/r;
    }
  	// Package conversions into a vector and Update KF
  	VectorXd zed(3);
  	zed << r, th, r_dot;
  	// Update KF
  	Update(zed);
}
