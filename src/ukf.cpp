#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initialization status
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // state dimension
  n_x_ = 5;

  // augmented state dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // time
  time_us_ = 100000;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .5; // need to modify

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .5; // need to modify

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

  // sigma point spreading parameter
  lambda_ = 3 - n_x_; // - n_aug_

  // set vector for weights_
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for(int i=1; i<2*n_aug_+1; i++)
  {
    weights_(i) = 0.5/(n_aug_+lambda_);
  }

  // NIS (Normalized Innovation Squared)
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if(!is_initialized_)
  {
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_.coeff(0);
      double phi = meas_package.raw_measurements_.coeff(1);
      double rho_dot = meas_package.raw_measurements_.coeff(2);

      x_ << rho*cos(phi), rho*sin(phi), rho_dot, 0, 0; // initialize x (radar)
    } else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << meas_package.raw_measurements_.coeff(0),
            meas_package.raw_measurements_.coeff(1), 0, 0, 0;
    }

    P_ <<
        1, 0, 0, 0, 0, // px
        0, 1, 0, 0, 0, // py
        0, 0, 10, 0, 0, // v
        0, 0, 0, 10, 0, // yaw
        0, 0, 0, 0, 0; // yaw_rate

    time_us_ = meas_package.timestamp_; // initialize timestamp
    is_initialized_ = true;
    return;
  }

  // Prediction
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  //cout << "predicted v = " << x_(2) << endl;

  // update
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR) UpdateRadar(meas_package);
  else UpdateLidar(meas_package);
  //cout << "updated v = " << x_(2) << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  VectorXd x_aug = VectorXd(n_aug_); // create augmented state
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_); // create augmented covariance
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1); // create sigma point matrix

  // 1. create augmented sigma points
  // augment x
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // augment P
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create augmented sigma points
  MatrixXd L = P_aug.llt().matrixL();
  Xsig_aug.col(0) = x_aug;
  for(int i=0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)*L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*L.col(i);
  }

  // 2. predict sigma points
  for(int i=0; i<2*n_aug_+1; i++)
  {
    // extract values
    double px_ = Xsig_aug(0,i);
    double py_ = Xsig_aug(1,i);
    double v_ = Xsig_aug(2,i);
    double yaw_ = Xsig_aug(3,i);
    double yawdot_ = Xsig_aug(4,i);
    double nu_a_ = Xsig_aug(5,i);
    double nu_yawdd_ = Xsig_aug(6,i);

    // predicted states
    double px_est, py_est;

    // avoid division by zero
    if(fabs(yawdot_) > 0.001)
    {
      px_est = px_ + v_/yawdot_*(sin(yaw_+yawdot_*delta_t) - sin(yaw_));
      py_est = py_ + v_/yawdot_*(cos(yaw_) - cos(yaw_+yawdot_*delta_t));
    }
    else
    {
      px_est = px_ + v_*delta_t*cos(yaw_);
      py_est = py_ + v_*delta_t*sin(yaw_);
    }

    double v_est = v_; // constant velocity
    //if(fabs(v_) > 20) cout << "v = " << v_ << endl;
    double yaw_est = yaw_ + yawdot_ * delta_t;
    double yawdot_est = yawdot_; // constant turn rate

    // add noise
    px_est += 0.5 * nu_a_ * delta_t * delta_t * cos(yaw_);
    py_est += 0.5 * nu_a_ * delta_t * delta_t * sin(yaw_);
    v_est += nu_a_ * delta_t;

    //cout << "nu_a_ = " << nu_a_ << endl;

    yaw_est += 0.5 * nu_yawdd_ * delta_t * delta_t;
    yawdot_est += nu_yawdd_ * delta_t;

    // write predicted sigma point
    Xsig_pred_(0, i) = px_est;
    Xsig_pred_(1, i) = py_est;
    Xsig_pred_(2, i) = v_est;
    Xsig_pred_(3, i) = yaw_est;
    Xsig_pred_(4, i) = yawdot_est;
  }

  // 3. predict state and covariance
  // predict state mean
  x_.fill(0);
  for(int i=0; i<2*n_aug_+1;i++)
  {
    x_ = x_ + Xsig_pred_.col(i) * weights_(i);
  }
  //cout << "x = " << x_ << endl;

  // predict state covariance matrix
  P_.fill(0);
  for(int i=0; i<2*n_aug_+1;i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3)); // angle normalization
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // transform sigma points to measurement space
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  for(int i=0; i<2*n_aug_+1; i++)
  {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  // predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i=0; i<2*n_aug_+1;i++)
  {
    z_pred += weights_(i) * Zsig.col(i);
  }
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S = S + R;

  VectorXd z = meas_package.raw_measurements_;
  VectorXd y = z - z_pred;

  // calculate laser NIS
  NIS_laser_ = y.transpose() * S.inverse() * y;

  // calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(this->n_x_, n_z);
  Tc.fill(0.0);
  for(int i=0; i<2*this->n_aug_ + 1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //NormalizeAngle(z_diff(1));
    VectorXd x_diff = this->Xsig_pred_.col(i) - this->x_;
    NormalizeAngle(x_diff(3));
    Tc += this->weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  this->x_ += K * y;
  this->P_ -= K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  // transform sigma points to measurement space
  for(int i=0; i<2*n_aug_+1; i++)
  {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig(0, i) = sqrt(px*px + py*py); // rho
    Zsig(1, i) = atan2(py, px); // phi
    Zsig(2, i) = (px*v1 + py*v2) / sqrt(px*px + py*py); // rhodot
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++)
  {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S += R;

  VectorXd z = meas_package.raw_measurements_;
  VectorXd y = z - z_pred;

  NIS_radar_ = y.transpose() * S.inverse() * y; // derive radar NIS

  // derive cross correlation matrix
  MatrixXd Tc = MatrixXd(this->n_x_, n_z);
  Tc.fill(0.0);
  for(int i=0; i<2*this->n_aug_+1;i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));

    VectorXd x_diff = this->Xsig_pred_.col(i) - this->x_; // state difference
    NormalizeAngle(x_diff(3));
    Tc += this->weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse(); // kalman gain
  VectorXd z_diff = z - z_pred;
  NormalizeAngle(z_diff(1));

  // update state mean and covariance matrix
  this->x_ += K*z_diff;
  this->P_ -= K*S*K.transpose();
}

void UKF::NormalizeAngle(double& phi){
  phi = atan2(sin(phi), cos(phi));
}
