#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.empty() || estimations.size() != ground_truth.size()) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;

  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    // calculate the error
    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication: calculate the squared error
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean squared error
  rmse /= estimations.size();

  //calculate the root mean squared error
  rmse = rmse.cwiseSqrt();

  //return the result
  return rmse;
}