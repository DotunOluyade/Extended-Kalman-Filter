#include "tools.h"
#include <iostream>
#include "kalman_filter.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}



void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_lidar_in, MatrixXd &R_radar_in, MatrixXd &Q_in) {
  x_ = x_in; // state vector
  P_ = P_in;  
  F_ = F_in; // State transition matrix
  H_ = H_in;  
  R_lidar = R_lidar_in;
  R_radar = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
     * KF Measurement update step
     */
    VectorXd y = z - H_ * x_;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_lidar;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    
    // new state
    x_ = x_ + (K * y);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
   /**
   * update the state by using Extended Kalman Filter equations
   */
 
    /**
     * KF Measurement update step
     */
    Tools tools;
    MatrixXd Hj = tools.CalculateJacobian(x_);
    /* Cartesian coordinates */
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    /*Cartesian to polar coordinate conversion steps */
    double rho = sqrt(px*px + py*py);
    double theta = atan2(py, px);
    double rho_dot = (px*vx + py*vy) / rho;
    VectorXd h = VectorXd(3);

    /* y = z - h(x)  function h(x) maps values from cartesian to polar coordinates 
     * function h(x) is the value h    
     */
    h << rho, theta, rho_dot; 
     
    VectorXd y = z - h; 
  
    /* the resulting angle phi in the y vector should be adjusted so that it is between -pi and pi*/
  	while ( y(1) > M_PI || y(1) < -M_PI ) {
      if ( y(1) > M_PI ) {
        std::cout<<"y(1) is greater than pi" << y(1);
        y(1) -= 2 * M_PI;  // substract 2Pi
      } else {
        std::cout<<"y(1) is less than pi" << y(1);
        y(1) += 2 * M_PI;  // add 2Pi
      }
  	}
  
    MatrixXd Ht = Hj.transpose();
    MatrixXd S = Hj * P_ * Ht + R_radar;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    // new state
    x_ = x_ + (K * y);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj) * P_;
  
}
