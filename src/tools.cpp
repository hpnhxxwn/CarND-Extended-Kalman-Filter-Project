#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    for (size_t i = 0; i < estimations.size(); i++) {
        VectorXd r = estimations[i] - ground_truth[i];
        r = r.array() * r.array();
        rmse += r;
    }

    rmse /= estimations.size();

    rmse = sqrt(rmse.array());

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
    MatrixXd Hj(3, 4);
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    if (fabs(px) < 0.0001 and fabs(py) < 0.0001){
        px = 0.0001;
        py = 0.0001;
    }
    float eq1 = px * px + py * py;

    if(fabs(eq1) < 0.0000001){
        eq1 = 0.0000001;
    }

    float eq2 = sqrt(eq1);
    float eq3 = eq1 * eq2;

    Hj << px / eq2, py / eq2, 0, 0,
         -py / eq1, px / eq1, 0, 0,
            py * (vx * py - vy * px) / eq3, px * (vy * px - vx * py) / eq3, px / eq2, py / eq2;

    return Hj;
};
