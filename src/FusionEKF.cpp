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
  Hj_      = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  
  // Measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  // Measurement matrix for radar is calculated dynamically
  
  noise_ax = 9;
  noise_ay = 9;
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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    
    // first measurement
    cout << "EKF Initialization: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd initial_x = Tools::PolarToCartesian(measurement_pack.raw_measurements_);
      
      // Since we're measuring with radar,
      // variance in velocity and distance are both going to be relatively low
      MatrixXd initial_P = MatrixXd(4,4);
      initial_P << 1,0,0,0,
                   0,1,0,0,
                   0,0,1,0,
                   0,0,0,1;

      // Set the transition matrix with dt = 0
      MatrixXd initial_F = Tools::TransitionMatrix(0);
      
      // Our initial H is going to be the jacobian matrix given our initial x
      MatrixXd initial_H = Tools::CalculateJacobian(initial_x);
      
      // Set initial measurement covariance to radar's values
      MatrixXd initial_R = R_radar_;
      
      // Set initial process covariance using given noise and dt = 0
      MatrixXd initial_Q = Tools::ProcessCovarianceMatrix(0, noise_ax, noise_ay);
      
      ekf_.Init(initial_x, initial_P, initial_F, initial_H, initial_R, initial_Q);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      float vx = 0;
      float vy = 0;
      
      // Our initial x is just the measurement with velocity 0.
      VectorXd initial_x = VectorXd(4);
      initial_x << px, py, vx, vy;
      
      // Initial P assumes a high uncertainty in velocity because we have not directly measured it.
      MatrixXd initial_P = MatrixXd(4,4);
      initial_P << 1,0,0,0,
                   0,1,0,0,
                   0,0,1000,0,
                   0,0,0,1000;
      
      // Set transition matrix with dt = 0;
      MatrixXd initial_F = Tools::TransitionMatrix(0);
      
      // Set initial H to the jacobian of our initial X
      MatrixXd initial_H = Tools::CalculateJacobian(initial_x);
      
      // Set initial measurement covariance to radar's values
      MatrixXd initial_R = R_laser_;
      
      // Set initial process covariance using given noise and dt=0
      MatrixXd initial_Q = Tools::ProcessCovarianceMatrix(0, noise_ax, noise_ay);
      
      ekf_.Init(initial_x, initial_P, initial_F, initial_H, initial_R, initial_Q);
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  //compute the time elapsed between the current and previous measurements
  // time difference expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  
  // Update state transition matrix with new dt
  ekf_.F_ = Tools::TransitionMatrix(dt);
  
  // Update process covariance matrix with time difference
  ekf_.Q_ = Tools::ProcessCovarianceMatrix(dt, noise_ax, noise_ay);
  
  //std::cout << "DT:" << endl << dt << endl;
  // Make prediction, but only if the time difference is enough
  if ( dt > 0.001 ) {
    previous_timestamp_ = measurement_pack.timestamp_;
    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.H_ = Tools::CalculateJacobian(ekf_.x_);
    
    // Extract measurement from the package
    VectorXd z = measurement_pack.raw_measurements_;
    
    ekf_.UpdateEKF(z);
    
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    
    // Extract measurement from the package
    VectorXd z = measurement_pack.raw_measurements_;
    
    ekf_.Update(z);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
