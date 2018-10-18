#include <iostream>
#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

// Calculate the RMSE from estimations and ground truths
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) {
    /*	
    	Pseudocode for R.M.S.E. Equation:

          RMSE = sqrt( 1/n sigma(1 to n of (xt(est) - xt(tru))^2 ) )

       	Conditions:

          1. Neither the estimation nor ground_truth vectors can have zero elements
          2. Both vectors must be the same shape

       	Possible Improvements:

          1. Make the size of the RMSE vector dependent on the vectors found within 
              the input vectors (currently hardcoded to 4)
          2. Change the way in which the estimations are being looped over (boost::combine and iterate)
    */
  	// Initiallize rmse to all zeros
  	VectorXd rmse(4);
  	rmse << 0, 0, 0, 0;
  
  	// Test for Proper Input Vector Sizes
	if(estimations.size() == 0 || ground_truth.size() == 0) {
    	std::cout << "Error - Estimations and Ground Truth Vectors must not be empty";
    	return rmse;
    }
    else if(estimations.size() != ground_truth.size()) {
    	std::cout << "Error - Estimations and Ground Truth Vectors must be the same length";
    	return rmse;
    }
  	// Create Sum of the Differences and Square this Sum
  	for(unsigned int i = 0; i < estimations.size(); ++i) {
      	// Difference Between Estimation and True Values
    	VectorXd diff = estimations[i] - ground_truth[i];
      	// Element-wise multiplication Squares each Difference
      	VectorXd diff_2 = diff.array() * diff.array();
      	rmse += diff_2;
    }
  	// Divide rmse by number of diff's computed
  	rmse /= estimations.size();
  	// Take the Square Root and return the result
  	return rmse.array().sqrt();
}

// Calculates the Jacobian from
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  	MatrixXd Hj(3, 4);
  	// Unpack elments of x_state
  	double px = x_state(0), 
  		   py = x_state(1), 
  		   vx = x_state(2), 
  		   vy = x_state(3);
  	// Perform and store frequent calucations first
  	double pSquared = px*px + py*py;
  	// if this is zero or very close to zero, initialize and return Hj
  	if(fabs(pSquared) < 0.0001) {
      	// Error Message
    	std::cout << "Error - Cannot Calculate Jacobian if px and py are both zero";
     	// Initialize Hj to all zeros
      	Hj << 0, 0, 0, 0,
      		0, 0, 0, 0,
      		0, 0, 0, 0;
      	return Hj; 
    }
  	double pSqrt = sqrt(pSquared);
  	double pxOverPSqrt = px / pSqrt;
  	double pyOverPSqrt = py / pSqrt;
  	double pSqrtCubed = pSqrt*pSqrt*pSqrt;
  
  	// Fill the Jacobian
  	Hj << pxOverPSqrt, pyOverPSqrt, 0, 0,
  		(-py/pSquared), px/pSquared, 0, 0,
  		py*(vx*py-vy*px)/pSqrtCubed, px*(vy*px-vx*py)/pSqrtCubed, pxOverPSqrt, pyOverPSqrt;
  
  	return Hj;
}
