#include <iostream>
#include <eigen3/Eigen/Eigen>

class RotatedPlaneManager
{
public:
   struct DataPoint
   {
      double theta;
      Eigen::Vector4d plane;
   };
public:
   void solveFromDataPoints(std::vector<DataPoint>& data_points);
   void computePlane(double theta, Eigen::Vector4d& plane);
protected:
   Eigen::Matrix<double, 4, 2> _A;
   Eigen::Matrix<double, 4, 1> _B;
};

void RotatedPlaneManager::solveFromDataPoints(std::vector<DataPoint>& data_points)
{
   const int N = 12 + data_points.size();

   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
   matrix.resize(
      4 * data_points.size(),
      12 + data_points.size());
   matrix.setZero();

   for(int i=0; i<data_points.size(); i++)
   {
      const double cos_theta = cos(data_points[i].theta);
      const double sin_theta = sin(data_points[i].theta);

      matrix(4*i+0, 0) = cos_theta;
      matrix(4*i+0, 4) = sin_theta;
      matrix(4*i+0, 8) = 1.0;
      matrix(4*i+0, 12+i) = -data_points[i].plane(0);

      matrix(4*i+1, 1) = cos_theta;
      matrix(4*i+1, 5) = sin_theta;
      matrix(4*i+1, 9) = 1.0;
      matrix(4*i+1, 12+i) = -data_points[i].plane(1);

      matrix(4*i+2, 2) = cos_theta;
      matrix(4*i+2, 6) = sin_theta;
      matrix(4*i+2, 10) = 1.0;
      matrix(4*i+2, 12+i) = -data_points[i].plane(2);

      matrix(4*i+3, 3) = cos_theta;
      matrix(4*i+3, 7) = sin_theta;
      matrix(4*i+3, 11) = 1.0;
      matrix(4*i+3, 12+i) = -data_points[i].plane(3);
   }

   auto svd = matrix.jacobiSvd(Eigen::ComputeThinV);
   Eigen::VectorXd solution = svd.matrixV().rightCols<1>();

   _A.leftCols<1>() = solution.segment<4>(0);
   _A.rightCols<1>() = solution.segment<4>(4);
   _B = solution.segment<4>(8);
}

void RotatedPlaneManager::computePlane(double theta, Eigen::Vector4d& plane)
{
   Eigen::Vector2d X;
   X << cos(theta), sin(theta);
   plane = _A*X + _B;
}

int main(int num_args, char** args)
{
   Eigen::Matrix3d R0;
   R0 <<
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0;

   Eigen::Vector3d O;
   O << 0.0, 0.0, 0.0;

   Eigen::Vector3d P;
   P << 0.0, 0.0, 0.0;

   Eigen::Vector3d n;
   n << 1.0, 0.0, 0.0;

   std::vector<RotatedPlaneManager::DataPoint> keypoints;
   for(int i=0; i<6; i++)
   {
      const double theta = double(i)*M_PI/6.0;

      Eigen::Matrix3d Rz;
      Rz <<
         cos(theta), -sin(theta), 0.0,
         sin(theta), cos(theta), 0.0,
         0.0, 0.0, 1.0;

      Eigen::Vector4d plane;
      plane.head<3>() = R0 * Rz * n;
      plane.tail<1>() = -( P.transpose() + O.transpose()*R0*Rz ) * n;

      RotatedPlaneManager::DataPoint pt;
      pt.theta = theta;
      pt.plane = plane;

      keypoints.push_back(pt);
   }

   RotatedPlaneManager manager;
   manager.solveFromDataPoints(keypoints);

   for(RotatedPlaneManager::DataPoint& pt : keypoints)
   {
      Eigen::Vector4d plane;
      manager.computePlane(pt.theta, plane);

      plane /= plane.head<3>().norm();
      pt.plane /= pt.plane.head<3>().norm();

      std::cout << plane.transpose() << std::endl;
      std::cout << pt.plane.transpose() << std::endl;
      std::cout << std::endl;
   }

   return 0;
}
