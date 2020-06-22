#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include<stdlib.h>
#include<Eigen/Core>
#include<Eigen/Dense>
#include"config.h"
#include"imuInitialization.h"
#include"rotation.h"
#include<cmath>

using namespace std;

int main(int argc, char** argv)
{
   clock_t time_stt=clock();
   if(argc!=2)
     {
       cout<<"usage: parameter_file"<<endl;
       return 1;
     }
/*
    Config::setParameterFile ( argv[1] ); //调用设置参数的函数
    //调用参数文件中的参数
    string dataset_dir = Config::get<string> ( "dataset_dir" );//返回Config->file_[dataset_dir]，也就是返回数据的绝对路径
    string output_dir = Config::get<string> ( "output_dir" );//返回Config->file_[doutput_dir]，结果保存的路径
    double w_ie= Config::get<double> ( "we" );
    double location[3]={Config::get<double> ( "longitude" ), Config::get<double> ( "latitude" ), Config::get<double> ( "height" )};
    //double g_parameter[6]={Config::get<double> ( "a1" ),Config::get<double> ( "a2" ),Config::get<double> ( "a3" ),Config::get<double> ( "a4" ),Config::get<double> ( "a5" ),Config::get<double> ( "a6" )};
    //读取观测数据
*/
    //调用参数文件中的参数
    cv::FileStorage fsSettings(argv[1], cv::FileStorage::READ);
    if(!fsSettings.isOpened())
    {
        std::cerr << "ERROR: Wrong path to settings" << std::endl;
    }
    string dataset_dir = fsSettings["dataset_dir"];//返回Config->file_[dataset_dir]，也就是返回数据的绝对路径
    string output_dir = fsSettings["output_dir"];//返回Config->file_[doutput_dir]，结果保存的路径
    double w_ie= fsSettings["we"];
    double location[3]={fsSettings["lonitude"], fsSettings["latitude"], fsSettings["height"]};
    //读取观测数据
    ifstream fin(dataset_dir);
    if ( !fin )
    {
        cout<<"please generate the associate file called imu.txt!"<<endl;
        return 1;
    }


// 1,第一步，进行初始化.............................

   //vector<double>acc_x,acc_y,acc_z,gyro_x,gyro_y,gyro_z;
   double sum_gyro_x = 0;
	double sum_gyro_y = 0;
	double sum_gyro_z = 0;
	double sum_acc_x = 0;
	double sum_acc_y = 0;	
	double sum_acc_z = 0;
   int count = 0; //用于控制读取的初始化数据
   //读取初始化部分数据
   while(!fin.eof()) //读到末尾前
      {
         string time,gyro_xs,gyro_ys,gyro_zs,acc_xs,acc_ys,acc_zs;
         double acc_x,acc_y,acc_z,gyro_x,gyro_y,gyro_z;         
         fin>>time>>gyro_xs>>gyro_ys>>gyro_zs>>acc_xs>>acc_ys>>acc_zs;
         gyro_x = atof(gyro_xs.c_str()); //c_str()是返回一个内容同sting类相同的C风格字符串
         gyro_y = atof(+gyro_ys.c_str());
         gyro_z = atof(gyro_zs.c_str());
         acc_x = atof(acc_xs.c_str());
	      acc_y = atof(acc_ys.c_str());
	      acc_z = atof(acc_zs.c_str());

	      sum_gyro_x += gyro_x;
	      sum_gyro_y += gyro_y;
	      sum_gyro_z += gyro_z;
         sum_acc_x  += acc_x;
	      sum_acc_y  += acc_y;
	      sum_acc_z  += acc_z;
         count++;
         if(count >=4000) break;//读取前4000个数据
       }
   fin.close();

   //求平均值
	sum_gyro_x = sum_gyro_x/count*200;
	sum_gyro_y = sum_gyro_y/count*200;
	sum_gyro_z = sum_gyro_z/count*200;
   sum_acc_x = -sum_acc_x/count*200;
	sum_acc_y = -sum_acc_y/count*200;
	sum_acc_z = -sum_acc_z/count*200;
     
	double acc_xyz[3] = {sum_acc_x,sum_acc_y,sum_acc_z};
	double gyro[3] = {sum_gyro_x,sum_gyro_y,sum_gyro_z};

   // 求惯导初始旋转矩阵
	Eigen::Matrix3d C_b_n=imuInitialization::getImuInitialization(w_ie,acc_xyz, gyro, location);
	cout<<C_b_n<<endl;
	double pich = -(atan(C_b_n(2,0)/sqrt(1-C_b_n(2,0)*C_b_n(2,0))))*180/M_PI;
	double roll = (atan2(C_b_n(2,1),C_b_n(2,2)))*180/M_PI;
	double yaw = (atan2(C_b_n(1,0),C_b_n(0,0)))*180/M_PI;
	cout<<"pich= "<<pich<<endl<<"roll= "<<roll<<endl<<"yaw= "<<yaw<<endl;
	cout<<"time use in QR compositiong is "<<1000*(clock()-time_stt)/(double) CLOCKS_PER_SEC<<"ms"<<endl;
	
	Eigen::Vector3d euler_angles = C_b_n.eulerAngles(2,1,0);
	cout<<"yaw pitch roll = "<<euler_angles.transpose()*180/M_PI<<endl<<endl;	

   //初始化参数
   //配置基本地球参数	

   Eth eth;
   /*
   eth.we = Config::get<double> ("we");
   eth.a = Config::get<double> ("a"); 
   eth.e = Config::get<double> ("e");
   //eth.e2 = Config::get<double> ("e2");
   */


   eth.we = fsSettings["we"];
   eth.a = fsSettings["a"]; 
   eth.e = fsSettings["e"];
   //eth.e2 = Config::get<double> ("e2");
   fsSettings.release();
   eth.e2 = sqrt(6378160*6378160 + 6356775*6356775)/6356775;
   eth.w_ie_e<<0,0,eth.we;
   
   Ins ins;
   ins.eth = eth;
   ins.pos<<location[0],location[1],location[2];
   ins.C_b_n = C_b_n;
   ins.q_b_n = Eigen::Quaterniond(C_b_n);//////初始旋转四元数
   ins.dv_f_k_n<<0,0,0;  //比力速度增量
   ins.dv_g_k_n<<0,0,0;  //有害速度增量
   ins.v_n<<0,0,0; //初始速度
   double lamd = ins.pos(0); double fi = ins.pos(1); //上一时刻的经度纬度
   
   Eigen::Quaterniond q_n_e(cos(-M_PI/4-fi/2)*cos(lamd/2),-sin(-M_PI/4-fi/2)*sin(lamd/2),sin(-M_PI/4-fi/2)*cos(lamd/2),cos(-M_PI/4-fi/2)*sin(lamd/2));  //2.79 
   ins.q_n_e = q_n_e;

   //读取数据
   /*
   string dataset_dir_data = Config::get<string> ( "dataset_dir_data" ); //imu.txt文件路径
   cout<<dataset_dir_data<<endl;
   ifstream fin1(dataset_dir_data);
    if ( !fin1 )
    {
        cout<<"please generate the associate file called associate.txt!"<<endl;
        return 1;
    }
    */

//2，开始更新................................
   vector< vector<Eigen::Vector3d> > output;
   vector<double> imu_time;
   bool ium_update = 0;
   fin.open(dataset_dir);
   while(!fin.eof()) //读到末尾前
      {
         string time,gyro_xs,gyro_ys,gyro_zs,acc_xs,acc_ys,acc_zs;
         double cur_dt,cur_acc_x,cur_acc_y,cur_acc_z,cur_gyro_x,cur_gyro_y,cur_gyro_z; //当前时刻数据
         double prev_dt,prev_acc_x,prev_acc_y,prev_acc_z,prev_gyro_x,prev_gyro_y,prev_gyro_z; //上一时刻数据         

         fin>>time>>gyro_xs>>gyro_ys>>gyro_zs>>acc_xs>>acc_ys>>acc_zs;
         cur_dt = atof(time.c_str());
         cur_acc_x = atof(gyro_xs.c_str());
         cur_acc_y = atof(gyro_ys.c_str());
         cur_acc_z = atof(gyro_zs.c_str());
         cur_gyro_x = atof(acc_xs.c_str());
	      cur_gyro_y = atof(acc_ys.c_str());
	      cur_gyro_z = atof(acc_zs.c_str());
         
         if(ium_update){
            //开始计算
            ins.dt = cur_dt - prev_dt;
            ins.wb_k<< cur_gyro_x, cur_gyro_y, cur_gyro_z;//当前时刻的角度增量观测值
            ins.wb_k1<< prev_gyro_x, prev_gyro_y, prev_gyro_z; //上一时刻的角度增量观测值
            ins.fb_k<< cur_acc_x, cur_acc_y, cur_acc_z; //当前时刻的速度增量观测值
            ins.fb_k1<< prev_gyro_x, prev_acc_y, prev_acc_z; //上一时刻的速度增量观测值

            //Ins& ins1=ins;
            //角速度更新
            ins = imuInitialization::extrapolation(ins);
            //速度更新
            ins = imuInitialization::getVelocity(ins);
            //位置更新
            ins = imuInitialization::getLocation(ins);
            //姿态更新
            ins = imuInitialization::getPosture(ins);
            //Plastance plastance=imuInitialization::getPlastance(ins.eth,ins.pos,ins.v_n);
            //保存结果

         
            Eigen::Vector3d angle;
            angle<<(double)-(atan(ins.C_b_n(2,0)/sqrt(1-ins.C_b_n(2,0)*ins.C_b_n(2,0))))*180/M_PI,
	                (double)(atan2(ins.C_b_n(2,1),ins.C_b_n(2,2)))*180/M_PI;
	                (double)(atan2(ins.C_b_n(1,0),ins.C_b_n(0,0)))*180/M_PI;
         
            vector<Eigen::Vector3d> result;
            result.push_back(ins.v_n);
            result.push_back(ins.pos);
            result.push_back(angle);
            output.push_back(result);
            imu_time.push_back(prev_dt);
         }
          ium_update = true;
          prev_dt = cur_dt;
          prev_acc_x = cur_acc_x;
          prev_acc_y = cur_acc_y;
          prev_acc_z = cur_acc_z;
          prev_gyro_x = cur_gyro_x;
          prev_gyro_y = cur_gyro_y;
          prev_gyro_z = cur_gyro_z; 

         if(fin.good()==false) break;//读到末尾强制退出
       }
       //输出结果
       
       for(int i = 0; i < output.size(); i++){
          ofstream fout(output_dir, ios::app);
          fout.setf(ios::fixed, ios::floatfield);
          fout.precision(5);
          //fout<<imu_time[i]<<output[i][0](0)<<" "<<output[i][0](1)<<" "<<output[i][0](2)<<"   "output[i][1](0)<<" "output[i][1](1)<<" "output[i][1](2)<<"   "output[i][2](0)<<" "output[i][2](1)<<" "output[i][2](2)<<endl
          fout << imu_time[i] <<" "
               << output[i][0].x() << " " << output[i][0].y() << " " << output[i][0].z() << " "
               << output[i][1].x() << " " << output[i][1].y() << " " << output[i][1].z() << " "
               << output[i][2].x() << " " << output[i][2].y() << " " << output[i][2].z() << " "
               << endl;
          fout.close();
       }

      return 0;
}
