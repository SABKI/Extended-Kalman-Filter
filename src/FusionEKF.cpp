#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>

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
  
  // measurement covariance
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  
  // measurement matrix
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

      /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  
  // state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;
  
  // initialize measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
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
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    // the state vector x_
    ekf_.x_ = VectorXd(4);
    cout << "EKF: " << endl;
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      
      // to convert from Polar Coordinates (rho,phi) to Cartesian Coordinates (x,y) :
      // positions  ---> x = rho × cos( phi ) & y = rho × sin( phi )
      // velocities  ---> vx = dx/dt = rho_dot × cos( phi ) & vy = dy/dt = rho_dot × sin( phi )
      
      //  let's first read from the radar data : rho_measured, phi_measured, rhodot_measured
      float rho_measured = measurement_pack.raw_measurements_[0]; 
      float phi_measured = measurement_pack.raw_measurements_[1]; 
      float rhodot_measured = measurement_pack.raw_measurements_[2]; 
      
      // important point when calculating the error y with radar sensor data is to normalize the angle Phi while calculating the error y vector 		so this angle would be between −π and +π by either adding or subtracting -+2π from \phi.
      
      while ( phi_measured > M_PI || phi_measured < -M_PI ) 
      {
           if ( phi_measured > M_PI ) {  phi_measured -= 2.0 * M_PI;} 
           else {  phi_measured += 2.0 * M_PI; }
      }
      
      // Polar to Cartesian
      double px  = rho_measured * cos(phi_measured);
      double py  = rho_measured * sin(phi_measured);
      double vx = rhodot_measured * cos(phi_measured);
      double vy = rhodot_measured * sin(phi_measured);
      
      ekf_.x_ << px,
                 py,
                 vx,
                 vy;
                                                    
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
    {
      // TODO: Initialize state.
      
         ekf_.x_ << measurement_pack.raw_measurements_[0], 
                    measurement_pack.raw_measurements_[1], 
                    measurement_pack.raw_measurements_[2], 
                    measurement_pack.raw_measurements_[3];   
    }

    
    // check that x and y are not both zero (the state vector will be used to calculate the Jacobian for the EKF --> avoid dividing by 0)
    if (fabs(ekf_.x_(0)) < 0.0001 && fabs(ekf_.x_(1)) <  0.0001)
    {
      ekf_.x_(0) =  0.0001;
      ekf_.x_(1) =  0.0001;
    }
    
    // Store the first timestamp_
    previous_timestamp_ = measurement_pack.timestamp_; 
    
    // Initial state setup
    cout << "Initial state of ekf_.x_: " << ekf_.x_ << endl;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
   // let's compute the time elapsed between the current and previous measurements, dt is expressed in seconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  // Update the transition matrix F_ so that the time is integrated
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // set the acceleration noise components
  float noise_ax = 9;
  float noise_ay = 9;
  // set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

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
    // TODO: Radar updates
    
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);  // In the radar update step, the Jacobian matrix Hj is used to calculate S, K and P
    ekf_.R_ = R_radar_; // we call the radar measurement covariance matrix
    ekf_.UpdateEKF(measurement_pack.raw_measurements_); // call the EKF Update function

  } else {
    // TODO: Laser updates
    
    ekf_.H_ = H_laser_; // we use H of the laser because we are in the linear motion model
    ekf_.R_ = R_laser_; // we call the Lidar measurement covariance matrix
    ekf_.Update(measurement_pack.raw_measurements_); // call the KF Update function

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
