#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <algorithm>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // if true, logs NIS values
  log_nis_ = true;

  // if true, log x_ and P_ values
  log_xP_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  // number of sigma columns
  n_sigma_col_ = 2 * n_aug_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  //create vector for weights
  weights_ = VectorXd(n_sigma_col_);
  weights_(0) = lambda_/(lambda_+n_aug_);

  for (int i = 1; i < n_sigma_col_; i++) {
    weights_(i) = 1/(2*(lambda_ + n_aug_));
  }

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  Xsig_pred_ = MatrixXd(n_x_, n_sigma_col_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // std_a_ = 2.0 && std_yawdd_ = 0.5 : 0.0711787
  // 0.0831147
  // 0.343067
  // 0.225984
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.


  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
  time_us_ = 0;

  // Set Q matrix
  Q_ = MatrixXd(2,2);
  Q_ << std_a_*std_a_, 0,
        0, std_yawdd_*std_yawdd_;

  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_*std_laspx_, 0,
             0, std_laspy_*std_laspy_;

  H_laser_ = MatrixXd(2,5);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;
  R_radar_ = MatrixXd(n_z_, n_z_);
  R_radar_.fill(0.0);
  R_radar_(0,0) = std_radr_*std_radr_;
  R_radar_(1,1) = std_radphi_*std_radphi_;
  R_radar_(2,2) = std_radrd_*std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  // initial measurement
  if(!is_initialized_) {
    /**
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
     **/
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
       * rho: meas_package.raw_measurements_[0]
       * phi: meas_package.raw_measurements_[1]
       * rho_dot: meas_package.raw_measurements_[2]
       * state -> x_ = px, py, v, psi, psi_dot
      */
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
//      double rho_dot = meas_package.raw_measurements_[2];

      x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
       * px: meas_package.raw_measurements_[0]
       * py: meas_package.raw_measurements_[1]
       * state -> x_ = px, py, v, psi, psi_dot
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // done initializing, no need to predict or update
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // prediction
  // compute the time elapsed between the current and previous measurements
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  UKF::Prediction(dt);

  // update
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UKF::UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UKF::UpdateLidar(meas_package);
  }

  if(log_xP_) {
    // print the output
    cout << "x_ = " << x_ << endl;
    cout << "P_ = " << P_ << endl;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  ///* Augmentation Step
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << x_, 0, 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_col_);

  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  ///* Prediction Step
  //predict sigma points
  for(int i = 0; i < n_sigma_col_; i++) {
    // capture augmented state values for each column
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double psi = Xsig_aug(3, i);
    double psi_dot = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_psidd = Xsig_aug(6, i);

    // calculate prediction state
    double px_p, py_p;

    //avoid division by zero
    if(fabs(psi_dot) < 0.00001) {
      px_p = px + v * cos(psi) * delta_t;
      py_p = py + v * sin(psi) * delta_t;
    } else {
      px_p = px + (v/psi_dot) * (sin(psi + psi_dot * delta_t) - sin(psi));
      py_p = py + (v/psi_dot) * (cos(psi) - cos(psi + psi_dot * delta_t));
    }

    double v_p = v;
    double psi_p = psi + psi_dot * delta_t;
    double psi_dot_p = psi_dot;

    // add noise
    px_p = px_p + 0.5 * delta_t*delta_t * cos(psi) * nu_a;
    py_p = py_p + 0.5 * delta_t*delta_t * sin(psi) * nu_a;
    v_p = v_p + delta_t * nu_a;
    psi_p = psi_p + 0.5 * delta_t*delta_t * nu_psidd;
    psi_dot_p = psi_dot_p + delta_t * nu_psidd;

    //write predicted sigma points into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = psi_p;
    Xsig_pred_(4, i) = psi_dot_p;
  }

  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sigma_col_; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sigma_col_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization. Make sure to
    // always normalize when calculating the
    // difference between angles
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  if(log_nis_) {
    // Calculate Normalized Innovation Squared (NIS)
    laser_nis_.push_back(y.transpose() * Si * y);
    // laser has 2 degrees of freedom
    double n_laser_over = count_if(laser_nis_.begin(), laser_nis_.end(), [](int i){ return i>5.991;});
    cout << "total laser measurements: " << laser_nis_.size()
         << " number over chi-squared value: " << n_laser_over
        << ", ratio: " << n_laser_over/laser_nis_.size() << endl;
  }

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, n_sigma_col_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);

  //transform sigma points into measurement space
  for(int i = 0; i < n_sigma_col_; i++) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);
    double psi_dot = Xsig_pred_(4,i);

    double rho = sqrt(px*px + py*py);
    double phi = atan2(py, px);
    double rho_dot = (px*cos(psi)*v + py*sin(psi)*v)/rho;

    Zsig.col(i) << rho, phi, rho_dot;
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i = 0; i < n_sigma_col_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  S.fill(0.0);
  for(int i = 0; i < n_sigma_col_; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_radar_;

  ///* Update the state and covariance matrices
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  //calculate cross correlation matrix
  for(int i = 0; i < n_sigma_col_; i ++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  if(log_nis_) {
    //Calculate Normalized Innovation Squared (NIS)
    radar_nis_.push_back(z_diff.transpose() * S.inverse() * z_diff);
    // radar has 3 degrees of freedom
    double n_radar_over = count_if(radar_nis_.begin(), radar_nis_.end(), [](int i){ return i>7.815;});
    cout << "total radar measurements: " << radar_nis_.size()
         << " number over chi-squared value: " << n_radar_over
        << ", ratio: " << n_radar_over/radar_nis_.size() << endl;
  }

  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

}
