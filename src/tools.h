#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  static VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  static MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
  * Calculates the Q matrix from the given timestamp and noise values)
  **/
  static MatrixXd ProcessCovarianceMatrix(float dt, float noise_ax, float noise_ay);
    
  static MatrixXd TransitionMatrix(float time_difference);
    
  static VectorXd PolarToCartesian(const VectorXd &z);
    
  static VectorXd CartesianToPolar(const VectorXd &x);
    
};

#endif /* TOOLS_H_ */
