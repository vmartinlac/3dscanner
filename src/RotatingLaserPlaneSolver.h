#pragma once

#include <eigen3/Eigen/Eigen>
#include <vector>

class RotatingLaserPlaneSolver
{
public:
   struct Example
   {
      double theta;
      Eigen::Vector4d plane;
   };
public:
   void learnFromExamples(std::vector<Example>& examples);
   void computeLaserPlane(double theta, Eigen::Vector4d& plane);
protected:
   Eigen::Matrix<double, 4, 2> _A;
   Eigen::Matrix<double, 4, 1> _B;
};

