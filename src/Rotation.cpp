#include"rotation.h"

Eigen::Matrix3d Rotation:: vector2anti_matrix(const Eigen::Vector3d& r_k)
  {
    Eigen::Matrix3d a;
    a<<0,-r_k(1),r_k(2),
       r_k(1),0,-r_k(3),
       -r_k(2),r_k(3),0;
    return a;
  };//向量转反对称矩阵
Eigen::Vector4d Rotation::qmul(const Eigen::Vector4d& q1,const Eigen::Vector4d& q2) //四元数乘法
  {  
    Eigen::Vector4d q;
        q<<q1(0) * q2(0) - q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3),
          q1(0) * q2(1) + q1(1) * q2(0) + q1(2) * q2(3) - q1(3) * q2(2),
          q1(0) * q2(2) + q1(2) * q2(0) + q1(3) * q2(1) - q1(1) * q2(3),
          q1(0) * q2(3) + q1(3) * q2(0) + q1(1) * q2(2) - q1(2) * q2(1);
    return q;
  }

Eigen::Matrix3d Rotation::quaternion2rotationMatrix(const Eigen::Quaterniond& q) //四元数转旋转矩阵
  {
     Eigen::Matrix3d rotationMtrix;
    /* rotationMtrix<<q(0)*q(0)+q(1)*q(1)-q(2)*q(2)-q(3)*q(3),2*(q(1)*q(2)-q(0)*q(3)),2*(q(1)*q(3)+q(0)*q(2)),
                    2*(q(1)*q(2)+q(0)*q(3)),q(0)*q(0)-q(1)*q(1)+q(2)*q(2)-q(3)*q(3),2*(q(2)*q(3)-q(0)*q(1)),
                    2*(q(1)*q(3)-q(0)*q(2)),2*(q(2)*q(3)+q(0)*q(1)),q(0)*q(0)-q(1)*q(1)-q(2)*q(2)+q(3)*q(3);
     */
     rotationMtrix<<q.w()*q.w()+q.x()*q.x()-q.y()*q.y()-q.z()*q.z(),2*(q.x()*q.y()-q.w()*q.z()),2*(q.x()*q.z()+q.w()*q.y()),
                    2*(q.x()*q.y()+q.w()*q.z()),q.w()*q.w()-q.x()*q.x()+q.y()*q.y()-q.z()*q.z(),2*(q.y()*q.z()-q.w()*q.x()),
                    2*(q.x()*q.z()-q.w()*q.y()),2*(q.y()*q.z()+q.w()*q.x()),q.w()*q.w()-q.x()*q.x()-q.y()*q.y()+q.z()*q.z();
     return rotationMtrix;
  }

Eigen::Quaterniond Rotation::eulerAngle2quaternion(const Eigen::Vector3d& eulerAngle) //欧拉角转四元数
  {
    double q0=cos(eulerAngle(0)/2)*cos(eulerAngle(1)/2)*cos(eulerAngle(2)/2)+sin(eulerAngle(0)/2)*sin(eulerAngle(1)/2)*sin(eulerAngle(2)/2);
    double q1=sin(eulerAngle(0)/2)*cos(eulerAngle(1)/2)*cos(eulerAngle(2)/2)-cos(eulerAngle(0)/2)*sin(eulerAngle(1)/2)*sin(eulerAngle(2)/2);
    double q2=cos(eulerAngle(0)/2)*sin(eulerAngle(1)/2)*cos(eulerAngle(2)/2)+sin(eulerAngle(0)/2)*cos(eulerAngle(1)/2)*sin(eulerAngle(2)/2);
    double q3=cos(eulerAngle(0)/2)*cos(eulerAngle(1)/2)*sin(eulerAngle(2)/2)-sin(eulerAngle(0)/2)*sin(eulerAngle(1)/2)*cos(eulerAngle(2)/2);
    Eigen::Quaterniond q(q0,q1,q2,q3);
    return q;
  }
