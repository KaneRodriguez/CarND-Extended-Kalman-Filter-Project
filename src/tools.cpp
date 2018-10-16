#include <iostream>
#include "tools.h"

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
  	for(int i = 0; i < estimations.size(); ++i) {
      	// Difference Between Estimation and True Values
    	auto diff = estimations[i] - ground_truth[i];
      	// Element-wise multiplication Squares each Difference
      	rmse += diff.array() * diff.array();
    }
  	// Divide rmse by number of diff's computed
  	rmse /= estimations.size();
  	// Take the Square Root and return the result
  	return rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
}
