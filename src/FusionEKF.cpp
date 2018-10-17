#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.F_ = MatrixXd(4, 4);
	ekf_.Q_ = MatrixXd(4, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
          0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
          0, 0.0009, 0,
          0, 0, 0.09;

	H_laser_ << 1, 0, 0, 0,
  				0, 1, 0, 0;
  
  	// Initialize the Transition Matrix
  	ekf_.F_ << 1, 0, 1, 0,
		  	   0, 1, 0, 1,
		       0, 0, 1, 0,
		       0, 0, 0, 1;
  	// Initialize the State Covariance Matrix 
  	ekf_.F_ << 1, 0, 0, 0,
		  0, 1, 0, 0,
		  0, 0, 1000, 0,
		  0, 0, 0, 1000;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // Initialize the state ekf_.x_ with the first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    	// Unpackage Raw Measurements (Radar measures in polar coordinates)
      	float r, th, r_dot = measurement_pack.raw_measurements_(0), 
      						 measurement_pack.raw_measurements_(1), 
     	 					 measurement_pack.raw_measurements_(2);
      	// Precalculation
      	float cosTh, sinTh = cos(th), sin(th);
      	// Convert radar from polar to cartesian coordinates and initialize state.
      	ekf_.x_ << r*cosTh, r*sinTh, r_dot*cosTh, r_dot*sinTh;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      	// Assign the measured px and py only. Laser (LIDAR) does not take velocity measurements.
      	ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      	ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    }
	// Initialize previous time stamp
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
    
    // Calculate Change in Time Since Last Time Stamp (in seconds)
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    // Update Time Stamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // Account for Noise From Acceleration Components
   	float ax = 9;
    float ay = 9;
  
    // Precalculations
    float dt_2 = dt*dt;
    float dt_3 = dt*dt_2;
    float dt_4 = dt*dt_3;
	float dt_3_ax_over_2 = dt_3*ax/2;
	float dt_3_ay_over_2 = dt_3*ay/2;
  
    // Update the state transition matrix F according to the new elapsed time
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

 	// Update the process noise covariance matrix
    ekf_.Q_ << dt_4/4*ax, 0, dt_3_ax_over_2, 0,
               0, dt_4/4*ay, 0, dt_3_ay_over_2,
               dt_3_ax_over_2, 0, dt_2*ax, 0,
               0, dt_3_ay_over_2, 0, dt_2*ay;
    
  	ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

	// Perform the Update Step Based on Sensor Type
  	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    	// Radar updates for state and covariance matrices
      	H_ = Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
      	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
      	// Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_radar_;
      	ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
