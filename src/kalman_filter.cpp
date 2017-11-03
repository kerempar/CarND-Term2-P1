#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;  // object state
  P_ = P_in;  // object covariance matrix
  F_ = F_in;  // state transition matrix
  H_ = H_in;  // measurement matrix
  R_ = R_in;  // measurement covariance matrix
  Q_ = Q_in;  // process covariance matrix
}

// This is the prediction step.
// Make sure to set F and Q before doing this.
void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  
  x_ = F_ * x_;
  P_ = (F_ * P_ * F_.transpose()) + Q_;
}

// This is the update function for a lidar.
// Make sure to set P, H, and R before doing this.
void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  
  MatrixXd Ht = H_.transpose();
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  // New Estimates
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}

// This is the update function for a radar.
// Make sure to set P, H, and R before doing this.
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  MatrixXd Ht = H_.transpose();
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  
  VectorXd y = z - Tools::CartesianToPolar(x_);
  
  // normalizing angles
  while (y(1) > M_PI)
  {
    y(1) -= 2 * M_PI;
  }
  while (y(1) < -M_PI)
  {
    y(1) += 2 * M_PI;
  }

  
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  // New Estimates
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}
