#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	VectorXd RMSE(4);
	RMSE << 0,0,0,0;
	// Checking validity of estimation and ground truth data
	if (estimations.size() == 0 || estimations.size() != ground_truth.size())
	{
		cout << "Invalid estimation or ground data size " << endl ;
		return RMSE;
	}

	for (int i = 0 ; i < estimations.size();++i)
	{
		VectorXd error = estimations[i] - ground_truth[i];

		// vector element -wise multiplication
		error = error.array()*error.array();
		RMSE += error ;
	}

	// Mean Calculation
	RMSE = RMSE/estimations.size();

	// Root square
	RMSE = RMSE.array().sqrt();
	return RMSE;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3,4);

	float px,py,vx,vy ;
	px = x_state(0);
	py = x_state(1);
	vx = x_state(2);
	vy = x_state(3);

	float c1,c2,c3 ;
	c1 = px*px + py*py ;
	c2 = sqrt(c1);
	c3 = c1*c2;

	// Division by zero
	if (fabs(c1) < 0.0001){
		cout << " Error-Calculate Jacobian() denominator is zero " << endl;
		return Hj;
	}

	// Jacobian Matrix

	Hj << (px/c2), (py/c2),0,0,
			-(py/c1),(px/c1),0,0,
			py*(vx*py - vy*px)/c3, px*(px*vy -py*vx)/c3 , px/c2, py/c2 ;
	return Hj;


}
