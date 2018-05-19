#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // Set initialization to false to start.
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // TODO: Check to see that the conditional statements for using laser and radar are functioning correctly. Got an error when I set this to false.

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
    //TODO: Tweak this
  std_a_ = 5; // Definitely less than 9 for cars.

  // Process noise standard deviation yaw acceleration in rad/s^2
    //TODO: Tweak this
  std_yawdd_ = M_PI; // 2pi is too high, so starting by splitting the difference between 2*pi and 0 by just going with pi.
  
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
  
  /**
  Initialize
  */
  n_x_ = 5;

  // initial state vector
  x_ = Eigen::VectorXd(n_x_);

  // initial covariance matrix
  P_ = Eigen::MatrixXd(n_x_, n_x_);

  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  P_.setIdentity(); // Initialize P as an identity matrix.

  // NIS variables
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
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
  if (!is_initialized_) {
    // First measurement
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_) {
        // Initialize x_, P_, previous_time, and anything else needed.
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];

      // Convert from polar to cartesian
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      float v = sqrt(vx * vx + vy * vy);

      x_ << px, py, v, 0.5015, 0.3528; // Start with noise vals from EKF.

    } else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_) {
      // Zeros for velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }

  /***********
   * Run predict and update.
   ************/
  delta_t_ = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  Prediction(delta_t_);

  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_){
      UpdateLidar(meas_package);
  }
  else if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_){
      UpdateRadar(meas_package);
  }
  previous_timestamp_ = meas_package.timestamp_;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
   * Estimate the object's location. Modify the state
   * vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /***
   * Predict the sigma points
   */

  //create sigma point matrix
  Xsig_ = Eigen::MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  A_ = P_.llt().matrixL();

  //set lambda for non-augmented sigma points
  lambda_ = 3 - n_x_;

  //set first column of sigma point matrix
  Xsig_.col(0) = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig_.col(i + 1) = x_ + sqrt(lambda_ + n_x_) * A_.col(i);
    Xsig_.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A_.col(i);
  }

  /****
   * Augment sigma points predicted in last step
   */

  // Create augmented mean state
  x_aug_ = Eigen::VectorXd(n_aug_);
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  // Create augmented covariance matrix
  P_aug_ = Eigen::MatrixXd(7, 7);
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5, 5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  // Create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  // Create augmented sigma points
  Xsig_aug_ = Eigen::MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++) {
      Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
      Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  /******
   * Predict sigma points
   ******/
  for (int i = 0; i< 2*n_aug_+1; i++){
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * (sin (yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_ = Eigen::MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
    }


    /*****
     * Predict Mean and Covariance
     *****/
    // set weights
    weights_ = Eigen::VectorXd(2*n_aug_+1);
    double weight_0 = lambda_ / (lambda_ + n_aug_);

    weights_(0) = weight_0;

    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }

    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        // state difference
        VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
        while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

        P_ = P_ + weights_(i) * x_diff_ * x_diff_.transpose();
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

   Mapping from your state vector to observation space in lidar is linear, like EKF.
   No reason to use any non-linear techniques, can use regular Kalman Filter equations at beginning of EKF lessons.

   H matrix will look a bit different than lessons; has only 4 components in UKF due to CTRV model.

  */
  // Get measurements from measurement package
  VectorXd z_ = meas_package.raw_measurements_;

  // Dimensions of measurement - only x and y for lidar
  int n_z_ = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig_.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurement model
    Zsig_(0, i) = p_x;
    Zsig_(1, i) = p_y;
  }

  // Mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  // Measurement covariance matrix
  MatrixXd S_ = Eigen::MatrixXd(n_z_, n_z_);
  S_.fill(0.0);
  for(int i = 0; i < 2*n_aug_ + 1; i++) {
    // Calculate the residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
    // Set the S matrix
    S_ = S_ + weights_(i) * z_diff_ * z_diff_.transpose();
  }

  // Add measurement noise
  MatrixXd R_ = Eigen::MatrixXd(n_z_, n_z_);
  R_ << std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;
  S_ = S_ + R_;

  // Cross-correlation matrix
  MatrixXd Tc_ = Eigen::MatrixXd(n_x_, n_z_);
  Tc_.fill(0.0);

  // Do the update
  // Loop through Z sigma point values
  for(int i = 0; i < 2*n_aug_ + 1; i++) {
    // Compute the residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;

    // Compute the state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;

    // Update the cross-correlation matrix
    Tc_ = Tc_ + weights_(i) * x_diff_ * z_diff_.transpose();
  }

  // Kalman gain matrix
  MatrixXd K_ = Tc_ * S_.inverse();

  // Calculate z residual
  VectorXd z_diff_ = z_ - z_pred_;

  // Calculate NIS
  NIS_laser_ = z_diff_.transpose() * S_.inverse() * z_diff_;

  // Update the state mean and covariance matrix
  x_ = x_ + K_ * z_diff_;
  P_ = P_ - K_ * S_ * K_.transpose();
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

   Predict sigma points, predict unit covariance, then update
  */
  /***
   * Predict radar sigma points
   */
  // Get measurements from measurement package
  VectorXd z_ = meas_package.raw_measurements_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  // Matrix for Z sigma points
  MatrixXd Zsig_ = Eigen::MatrixXd(n_z_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred_ = Eigen::VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S_ = Eigen::MatrixXd(n_z_,n_z_);
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R_ = Eigen::MatrixXd(n_z_,n_z_);
  R_ <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S_ = S_ + R_;

  /*****
   * Update radar
   ****/
  //create example vector for incoming radar measurement
  z_ = VectorXd(n_z_);

  //create matrix for cross correlation Tc
  MatrixXd Tc_ = MatrixXd(n_x_, n_z_);
  Tc_.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    //residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
    //angle normalization
    while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

    Tc_ = Tc_ + weights_(i) * x_diff_ * z_diff_.transpose();
  }

  //Kalman gain K;
  MatrixXd K_ = Tc_ * S_.inverse();

  //residual
  VectorXd z_diff_ = z_ - z_pred_;

  //angle normalization
  while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
  while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

  // Calculate the NIS value
  NIS_radar_ = z_diff_.transpose() * S_.inverse() * z_diff_;

  //update state mean and covariance matrix
  x_ = x_ + K_ * z_diff_;
  P_ = P_ - K_ * S_ * K_.transpose();

}
