#include "kalman_filter.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;


// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {

	x_ = F_*x_;
	MatrixXd F_t;
	F_t = F_.transpose();
	// Evaluating State Co-variance Matrix
	P_ = F_*P_*F_t + Q_ ;

}

void KalmanFilter::Update(const VectorXd &z) {

	VectorXd z_pred = H_*x_;
	VectorXd y = z - z_pred ;
	MatrixXd H_t = H_.transpose();
	MatrixXd S = H_*P_*H_t + R_ ;
	MatrixXd S_inv = S.inverse();
	MatrixXd PH_t = P_*H_t;

	// Kalman Filter Gain
	MatrixXd K = PH_t*S_inv;

	// New Estimation of State and State Covariance Matrix
	x_ = x_ + (K*y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size,x_size);
	P_ = (I - K*H_)*P_;


}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

	float rho, phi,rho_dot;
	float px,py,vx,vy;
	px = x_[0];
	py = x_[1];
	vx = x_[2];
	vy = x_[3];
	rho = sqrt((px*px+py*py));
	phi = atan2(py,px);
	rho_dot = (px*vx+py*vy)/rho;


	VectorXd z_pred(3);
	z_pred << rho,phi,rho_dot;
	VectorXd y = z - z_pred ;
	// Making sure phi is between -pi and pi
	y(1) = atan2(sin(double(y(1))),cos(double(y(1))));



	MatrixXd H_t = H_.transpose();

	MatrixXd S = H_*P_*H_t + R_ ;

	MatrixXd S_inv = S.inverse();
	MatrixXd PH_t = P_*H_t;
	// Kalman Filter Gain
	MatrixXd K = PH_t*S_inv;
	// New Estimation of State and State Covariance Matrix
	x_ = x_ + (K*y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size,x_size);
	P_ = (I - K*H_)*P_;


}
