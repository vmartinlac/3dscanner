#include <iostream>
#include "RotatingLaserPlaneSolver.h"

void reference(double theta, Eigen::Vector4d& plane)
{
   Eigen::Matrix3d R0;
   R0 <<
      0.0, 1.0, 0.0,
      1.0, 0.0, 0.0,
      0.0, 0.0, 1.0;

   Eigen::Vector3d O;
   O << 0.0, 1.0, 2.0;

   Eigen::Vector3d P;
   P << -1.0, 2.0, 0.0;

   Eigen::Vector3d n;
   n << 1.0, 1.0, 0.0;

   Eigen::Matrix3d Rz;
   Rz <<
      cos(theta), -sin(theta), 0.0,
      sin(theta), cos(theta), 0.0,
      0.0, 0.0, 1.0;

   plane.head<3>() = R0 * Rz * n;
   plane.tail<1>() = -( P.transpose() + O.transpose()*R0*Rz ) * n;
}

int main(int num_args, char** args)
{
   std::vector<RotatingLaserPlaneSolver::Example> examples;
   for(int i=0; i<6; i++)
   {
      RotatingLaserPlaneSolver::Example pt;
      pt.theta = double(i)*M_PI/6.0;
      reference(pt.theta, pt.plane);

      examples.push_back(pt);
   }

   RotatingLaserPlaneSolver manager;
   manager.learnFromExamples(examples);

   for(int i=0; i<10; i++)
   {
      const double theta = double(i)*M_PI/3.0;

      Eigen::Vector4d plane;
      manager.computeLaserPlane(theta, plane);

      Eigen::Vector4d ref;
      reference(theta, ref);

      plane /= plane.head<3>().norm();
      ref /= ref.head<3>().norm();

      std::cout << plane.transpose() << std::endl;
      std::cout << ref.transpose() << std::endl;
      std::cout << std::endl;
   }

   return 0;
}
