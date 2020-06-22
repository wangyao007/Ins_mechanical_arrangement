#ifndef IMUINITIALIZATION_H
#define IMUINITIALIZATION_H

#include<iostream>
#include<string>
#include"config.h"
#include"rotation.h"
#include<Eigen/Core>
#include<Eigen/Dense>
#include<cmath>
#include <memory>


struct Eth //用于存放基本地球参数
{
 public:
  double e; //地球第一离心率
  double e2; // 地球第二离心率
  double a;  // 地球长半轴
  double we; // 地球自转角速度
  Eigen::Vector3d w_ie_e;
}; 

struct Plastance //用于存放角速度
{
  public:
   Eigen::Vector3d w_ie_n; 
   Eigen::Vector3d w_en_n;
   Eigen::Vector3d w_in_n;
   Eigen::Vector3d gn; //重力加速
};
struct Ins
{
  public:
   Eth eth;  //地球参数
   Eigen::Vector3d pos; //当前时刻位置
   Eigen::Matrix3d C_b_n; //当前旋转矩阵
   Eigen::Quaterniond q_b_n; //当前旋转四元数
   Eigen::Vector3d dv_f_k_n; //当前比力速度增量
   Eigen::Vector3d dv_g_k_n; //当前有害速度增量 
   Eigen::Vector3d v_n;  // 当前速度
   Eigen::Quaterniond q_n_e; //当前n系到e系的旋转四元数
   double dt;  //当前时间差
   Eigen::Vector3d wb_k; //当前陀螺测量增量
   Eigen::Vector3d fb_k; //当前加速度计测量增量
   Eigen::Vector3d fb_k1; //上一时刻加速度计测量增量
   Eigen::Vector3d wb_k1; //上一时刻陀螺测量增量
   Plastance plastance;   //地球角速度
   Eigen::Vector3d pos_k12; // 外推中间时刻位置
   Eigen::Vector3d v_n_k12; // 外推中间时刻速度
   Eigen::Vector3d v_n_k1;  //上一时刻速度
   Eigen::Vector3d pos_k1;  //上一时刻位置
   Eigen::Vector3d att;  
};


class imuInitialization
{
   private:
    imuInitialization(){}
   public:   
    static  Eigen::Matrix3d getImuInitialization(const double& w_ie,const double acc_xyz[],const double gyro[],const double location[]); 
                                                                                                //求初始姿态矩阵
    static double NormalGrvaity(const double& latitude, const double& height); //求重力加速度
    static Plastance getPlastance(Eth& eth,Eigen::Vector3d& location,Eigen::Vector3d& v);//求角速度
    static Ins extrapolation(Ins& ins);  //求位置速度外推
    static Ins getVelocity(Ins& ins); //速度更新
    static Ins getLocation(Ins& ins); //位置更新
    static Ins getPosture(Ins& ins); //姿态更新
   ~imuInitialization();
};




#endif
