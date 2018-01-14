#include "RotatingLaserPlaneSolver.h"

void RotatingLaserPlaneSolver::learnFromExamples(std::vector<Example>& examples)
{
   const int N = 12 + examples.size();

   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
   matrix.resize(
      4 * examples.size(),
      12 + examples.size());
   matrix.setZero();

   for(int i=0; i<examples.size(); i++)
   {
      const double cos_theta = cos(examples[i].theta);
      const double sin_theta = sin(examples[i].theta);

      matrix(4*i+0, 0) = cos_theta;
      matrix(4*i+0, 4) = sin_theta;
      matrix(4*i+0, 8) = 1.0;
      matrix(4*i+0, 12+i) = -examples[i].plane(0);

      matrix(4*i+1, 1) = cos_theta;
      matrix(4*i+1, 5) = sin_theta;
      matrix(4*i+1, 9) = 1.0;
      matrix(4*i+1, 12+i) = -examples[i].plane(1);

      matrix(4*i+2, 2) = cos_theta;
      matrix(4*i+2, 6) = sin_theta;
      matrix(4*i+2, 10) = 1.0;
      matrix(4*i+2, 12+i) = -examples[i].plane(2);

      matrix(4*i+3, 3) = cos_theta;
      matrix(4*i+3, 7) = sin_theta;
      matrix(4*i+3, 11) = 1.0;
      matrix(4*i+3, 12+i) = -examples[i].plane(3);
   }

   auto svd = matrix.jacobiSvd(Eigen::ComputeThinV);
   Eigen::VectorXd solution = svd.matrixV().rightCols<1>();

   _A.leftCols<1>() = solution.segment<4>(0);
   _A.rightCols<1>() = solution.segment<4>(4);
   _B = solution.segment<4>(8);
}

void RotatingLaserPlaneSolver::computeLaserPlane(double theta, Eigen::Vector4d& plane)
{
   Eigen::Vector2d X;
   X << cos(theta), sin(theta);
   plane = _A*X + _B;
}

