#ifndef ROTATION_H
#define ROTATION_H

#include<Eigen/Core>
#include<Eigen/Dense>
#include<iostream>
#include<cmath>
#include <memory>

class Rotation
{ 
  private:
  public:
  static Eigen::Matrix3d vector2anti_matrix(const Eigen::Vector3d& r_k);//向量转反对称旋转矩阵
  static Eigen::Vector4d qmul(const Eigen::Vector4d& q1,const Eigen::Vector4d& q2);//
  static Eigen::Matrix3d quaternion2rotationMatrix(const Eigen::Quaterniond& q);
  static Eigen::Quaterniond eulerAngle2quaternion(const Eigen::Vector3d& eulerAngle);
};
#endif
