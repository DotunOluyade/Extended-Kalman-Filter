#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
  
   H_laser_ << 1,0,0,0,
  			   0,1,0,0;
   
 
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "the first measurement is received; "
         << "measurement is performed by the "
         << ((measurement_pack.sensor_type_ == MeasurementPackage::RADAR) ? "radar " : "lidar ")
         << "sensor"
         << endl;
    
    VectorXd x_ = VectorXd(4);
    x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      // and initialize state.
      double rho = measurement_pack.raw_measurements_[0];       //range
      double phi = measurement_pack.raw_measurements_[1];       //bearing
      double rho_dot = measurement_pack.raw_measurements_[2];   // velocity of rho
      //Coordinates convertion from polar to cartesian
      double x = rho * cos(phi); 
      double y = rho * sin(phi);       
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      
      x_ << x,y,vx,vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
      x_ << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],0.0,0.0;
    }
    
    // state covariance matrix (initial velocity is unknown, hence the level of uncertainty is high)
    MatrixXd P(4, 4);
    P << 1, 0,    0,    0,
         0, 1,    0,    0,
         0, 0, 1000,    0,
         0, 0,    0, 1000;

    // state transition matrix (initially Δt is 0)
    MatrixXd F(4, 4);
    F << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;

    // process covariance matrix (initially Δt is 0, hence Q consists of 0's;
    // Eigen initializes matrices with 0's by default)
    MatrixXd Q(4, 4);

    ekf_.Init(x_, P, F, H_laser_, R_laser_, R_radar_, Q);

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // State transition matrix update
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // measurement noise
  double noise_ax = 9.0;
  double noise_ay = 9.0;

  // Update the process covariance matrix
  
  double dt_2 = dt * dt;    //dt^2
  double dt_3 = dt_2 * dt;  //dt^3
  double dt_4 = dt_3 * dt;  //dt^4
  double dt_4_4 = dt_4 / 4; //dt^4/4
  double dt_3_2 = dt_3 / 2; //dt^3/2
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
	         0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
	         dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
 	         0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;
  
  // Prediction Step
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
  	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
  	ekf_.Update(measurement_pack.raw_measurements_);
  }
  
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
