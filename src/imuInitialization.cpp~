#include"imuInitialization.h"

  //Eigen::Matrix3d getImuInitialization(const double& w_ie const double acc_xyz[], const double gyro[],const double location[], const double& g_parameter[] )// 参数w_ie为地球自传角速度，acc_xyz[3]为加速度记测量平均值，gyro[3]为陀螺测量平均值 
               //location[3]为经纬度、高程，g_parameter[4]为重力参数

   //计算旋转矩阵函数
   Eigen::Matrix3d imuInitialization::getImuInitialization(const double& w_ie,const double acc_xyz[], const double gyro[],const double location[] )
  {
    double g=NormalGrvaity(location[1],location[2]);//求重力加速度
    Eigen::Vector3d r_n,r_b, w_ie_n,w_ie_b;
    r_n<<0,0,g;
    r_b<<acc_xyz[0],acc_xyz[1],acc_xyz[2];
    w_ie_n<<w_ie*cos(location[1]*M_PI/180),0,-w_ie*sin(location[1]*M_PI/180);
    w_ie_b<<gyro[0],gyro[1],gyro[2];
    Eigen::Vector3d v_n=r_n.cross(w_ie_n);
    Eigen::Vector3d v_b=r_b.cross(w_ie_b);

    //Eigen::Vector3d v_n=CrossProduct(r_n,w_ie_n);
    //Eigen::Vector3d v_b=CrossProduct(r_b,w_ie_b);

    Eigen::Matrix3d lilun,celiang;
    lilun<<0,0,g,
           (double)w_ie_n[0],(double)w_ie_n[1],(double)w_ie_n[2],
           (double)v_n[0],(double)v_n[1],(double)v_n[2];
    celiang<<(double)r_b[0],(double)r_b[1],(double)r_b[2],
             (double)w_ie_b[0],(double)w_ie_b[1],(double)w_ie_b[2],
             (double)v_b[0],(double)v_b[1],(double)v_b[2];  
    Eigen::Matrix3d c_n_b=lilun.colPivHouseholderQr().solve(celiang);//采用QR分解的方式求逆矩阵
    //Eigen::Matrix3d c_n_b=lilun.inverse()*celiang;
    return c_n_b;

/*    lilun<<0,(double)w_ie_n[0],(double)v_n[0],
           0,(double)w_ie_n[1],(double)v_n[1],
           g,(double)w_ie_n[2],(double)v_n[2];
       
    celiang<<(double)r_b[0],(double)w_ie_b[0],(double)v_b[0],
             (double)r_b[1],(double)w_ie_b[1],(double)v_b[1],
	     (double)r_b[2],(double)w_ie_b[2],(double)v_b[2];  
    Eigen::Matrix3d c_n_b=lilun*(celiang.inverse());
*/    return c_n_b;
  }

   //求重力参数函数   
double imuInitialization::NormalGrvaity(const double& latitude, const double& height)
  { 
   double a1=9.7803267715;
   double a2=0.0052790414;
   double a3=0.0000232718;
   double a4=-0.000003087691089;
   double a5=0.000000004397731;
   double a6=0.00000000000721;
   double s2=sin(latitude*M_PI/180)*sin(latitude*M_PI/180);
   double s4=s2*s2;
   double ng=a1*(1+a2*s2+a3*s4)+(a4+a5*s2)*height+a6*height*height;
   //std::cout<<"ng is "<<ng<<std::endl;
   return ng;
  }

//求角速度,该函数是已知位置和速度，求地球角速度
Plastance imuInitialization::getPlastance(Eth& eth,Eigen::Vector3d& location,Eigen::Vector3d& v)//求角速度
  {
    double h=location(2);
    double e=eth.e;
    double a=eth.a;
    double we=eth.we;
    double M=a*(1-e*e)/pow(1-e*e*sin(location(1))*sin(location(1)),3/2); 
           //地球子午圈曲率半径M=a*(1-e*e)/(1-e*e*sin(location(1))*sin(location(1)))^(3/2); 
    double N=a/sqrt(1-e*e*sin(location(1))*sin(location(1)));  //地球卯酉圈半径
/*    Eigen::Vector3d w_en_n;w_en_n<<v(1)/(N+h),-v(0)/(M+h),-v(1)*tan(location(1))/(N+h); 
    Eigen::Vector3d w_ie_n;w_ie_n<<we*cos(location(1)),0,-we*sin(location(1));
    Eigen::Vector3d gn;gn<<0,0,NormalGrvaity(location(1),location(2));
    Eigen::Vector3d w_in_n=w_en_n+w_ie_n;
    Plastance plastance(w_en_n,w_ie_n,gn);
*/   Plastance plastance; 
    plastance.w_en_n<<v(1)/(N+h),-v(0)/(M+h),-v(1)*tan(location(1))/(N+h); //w_en_n,n系相对于e系的旋转角速度在n系的投影;
    plastance.w_ie_n<<we*cos(location(1)),0,-we*sin(location(1));
    plastance.gn<<0,0,NormalGrvaity(location(1),location(2)); //计算重力加速度
    plastance.w_in_n=plastance.w_en_n+plastance.w_ie_n;
    //创建应用以便返回引用
    return plastance;
  }


//求位置速度外推
Ins imuInitialization::extrapolation(Ins& ins)  
  {
 /*   Eth& eth=ins.eth; 
    Eigen::Vector3d& pos=ins.pos;
    Eigen::Vector3d& v_n=ins.v_n;
    Plastance& plastance=getPlastance(eth,pos,v_n); //求角速度
*/  
    ins.plastance=getPlastance(ins.eth,ins.pos,ins.v_n); //求角速度
    Eigen::Vector3d tao=ins.plastance.w_in_n*0.5*ins.dt/2;  
    Eigen::Vector3d epu=ins.eth.w_ie_e*0.5*ins.dt/2;  

    //求绝对值
    double m_tao=sqrt(tao(0)*tao(0)+tao(1)*tao(1)+tao(2)*tao(2));
    double m_epu=sqrt(epu(0)*epu(0)+epu(1)*epu(1)+epu(2)*epu(2));

    Eigen::Quaterniond q1(cos(m_epu),-sin(m_epu)/m_epu*epu(0),-sin(m_epu)/m_epu*epu(1),-sin(m_epu)/m_epu*epu(2));  // 2.74
    double lamd=ins.pos(0);double fi=ins.pos(1); //上一时刻的经度纬度
    Eigen::Quaterniond q2(cos(-M_PI/4-fi/2)*cos(lamd/2),-sin(-M_PI/4-fi/2)*sin(lamd/2),sin(-M_PI/4-fi/2)*cos(lamd/2),cos(-M_PI/4-fi/2)*sin(lamd/2)); //2.79 
    Eigen::Quaterniond q3(cos(m_tao),sin(m_tao)/m_tao*tao(0),sin(m_tao)/m_tao*tao(1),sin(m_tao)/m_tao*tao(2));// 2.73
    Eigen::Quaterniond q=q1*(q2*q3);   //2.72四元数更新 

    ins.pos_k12(0)=2*atan(q.z()/q.x());   //四元数反推经度q(3)/q(1),注意四元数的访问为q.w(),q.x(),q.y(),q.z();
    ins.pos_k12(1)=-0.5*M_PI-2*atan(-q.x()/q.z());//四元数反推纬度q(1)/q(3)
    ins.pos_k12(2)=ins.pos(2)-ins.v_n(2)*ins.dt*0.5; //2.75  外推大地高
    
    ins.v_n_k12=ins.v_n+0.5*(ins.dv_f_k_n+ins.dv_g_k_n); //7.26速度外推，获取中间时刻速度
    ins.plastance=getPlastance(ins.eth,ins.pos_k12,ins.v_n_k12);  //再次更新角速度 
    return ins;
  }


 //速度更新
Ins imuInitialization::getVelocity(Ins& ins) 
  {
   Eigen::Matrix3d I;
   I=Eigen::Matrix3d::Identity(3,3); //初始化为3*3的单位矩阵
   Eigen::Vector3d r_k=ins.plastance.w_in_n*ins.dt;

   Eigen::Vector3d dv_fk_bk1=ins.fb_k+0.5*ins.wb_k.cross(ins.fb_k)+1/12*(ins.wb_k1.cross(ins.fb_k)+ins.fb_k1.cross(ins.wb_k));// (2.68)
   ins.dv_f_k_n=(I-Rotation::vector2anti_matrix(0.5*r_k))*ins.C_b_n*dv_fk_bk1;// (2.69) 改一下0.5的位置%%%%%%%%%%%%%%
  
   //代入公式，求有害加速
   ins.dv_g_k_n=ins.plastance.gn*ins.dt-((2*ins.plastance.w_ie_n).cross(ins.plastance.w_en_n)).cross(ins.v_n_k12)*ins.dt; //(2.78)

   // 求更新速度 
   ins.v_n_k1=ins.v_n;  //记录上一时刻的速度，用于位置更新的外推
   ins.v_n=ins.v_n+ins.dv_f_k_n+ins.dv_g_k_n;  //(2.66)
   return ins;
  }

 //位置更新
 Ins imuInitialization::getLocation(Ins& ins)
  {
    Eigen::Vector3d v_n_k12=0.5*(ins.v_n+ins.v_n_k1); //2.84求中间速度
    ins.plastance=getPlastance(ins.eth,ins.pos_k12,v_n_k12); // %位置还用之前外推的中间位置
    Eigen::Vector3d tao=ins.plastance.w_in_n*ins.dt/2;  
    Eigen::Vector3d epu=ins.eth.w_ie_e*ins.dt/2;  
    
    //求绝对值
    double m_tao=sqrt(tao(0)*tao(0)+tao(1)*tao(1)+tao(2)*tao(2));
    double m_epu=sqrt(epu(0)*epu(0)+epu(1)*epu(1)+epu(2)*epu(2));
    Eigen::Quaterniond q1(cos(m_epu),-sin(m_epu)/m_epu*epu(0),sin(m_epu)/m_epu*epu(1),sin(m_epu)/m_epu*epu(2));  //2.74
    Eigen::Quaterniond q2=ins.q_n_e;  //%%%%%%%%%%%%%采用更新的q_n_e，而不是位置结果在计算，否则会损失精度
    Eigen::Quaterniond q3(cos(m_tao),sin(m_tao)/m_tao*tao(0),sin(m_tao)/m_tao*tao(1),sin(m_tao)/m_tao*tao(2));  //2.73
    Eigen::Quaterniond q=q1*(q2*q3);   //2.72四元数更新 
    
    //四元数规范化
    double n2 = q.w()*q.w()+q.x()*q.x()+q.y()*q.y()+q.z()*q.z();  // q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3);
    if (n2>1.000001 || n2<0.999999)
      {
        double nq = 1/sqrt(n2); 
        q=Eigen::Quaterniond(q.w()*nq,q.x()*nq,q.y()*nq,q.z()*nq);
      }

    //保留结果
    ins.q_n_e=q; //将当期更新的q_e_n保留，用作下一时刻的位置更新
    ins.pos_k1=ins.pos;//保留上一时刻的位置，用于姿态更新计算中间时刻位置
    ins.pos(1)=2*atan((double)q.z()/q.w());   //四元数反推经度ins.pos(1)=2*atan((double)q(3)/q(0)); 
    ins.pos(2)=-0.5*M_PI-2*atan((double)-q.y()/q.z()); //四元数反推纬度ins.pos(2)=-0.5*M_PI-2*atan((double)-q(2)/q(4));
    ins.pos(3)=ins.pos(3)-v_n_k12(3)*ins.dt;
    return ins;
  }

   //姿态更新
Ins imuInitialization::getPosture(Ins& ins) 
  {
    Eigen::Vector3d r1_k=ins.wb_k+(1/12)*(ins.wb_k1,ins.wb_k);  //（2.103）
    double abs_r1_k=sqrt(r1_k(0)*r1_k(0)+r1_k(1)*r1_k(1)+r1_k(2)*r1_k(2));  //求绝对值
    Eigen::Quaterniond q_bk_bk1(cos(0.5*abs_r1_k),sin(0.5*abs_r1_k)/abs_r1_k*r1_k(0),sin(0.5*abs_r1_k)/abs_r1_k*r1_k(1),sin(0.5*abs_r1_k)/abs_r1_k*r1_k(2));    // (2.100)

    //根据位置更新结果，计算w_in_n_k12
    Eigen::Vector3d p=(ins.pos+ins.pos_k1)*0.5;  //先计算中间位置
    ins.plastance=getPlastance(ins.eth,p,ins.v_n_k12);

    //计算q_nk1_nk
    Eigen::Vector3d r2_k=ins.plastance.w_in_n*ins.dt;
    double abs_r2_k=sqrt(r2_k(0)*r2_k(0)+r2_k(1)*r2_k(1)+r2_k(2)*r2_k(2));   //求绝对值
    Eigen::Quaterniond q_nk1_nk(cos(0.5*abs_r2_k),-sin(0.5*abs_r2_k)/abs_r2_k*r2_k(0),-sin(0.5*abs_r2_k)/abs_r2_k*r2_k(1),-sin(0.5*abs_r2_k)/abs_r2_k*r2_k(2));   // (2.104)

    //计算四元数更新q_bk_nk
    Eigen::Quaterniond q_bk_nk=q_nk1_nk*(ins.q_b_n*q_bk_bk1); //(2.99)
   
    //四元数的规范化
    double n2 = q_bk_nk.w()*q_bk_nk.w()+q_bk_nk.x()*q_bk_nk.x()+q_bk_nk.y()*q_bk_nk.y()+q_bk_nk.z()*q_bk_nk.z();
    //double n2 = q_bk_nk(0)*q_bk_nk(0)+q_bk_nk(1)*q_bk_nk(1)+q_bk_nk(2)*q_bk_nk(2)+q_bk_nk(3)*q_bk_nk(3);
    if (n2>1.000001 || n2<0.999999)
     {
       double nq = 1/sqrt(n2); 
       q_bk_nk=Eigen::Quaterniond(nq*q_bk_nk.w(),nq*q_bk_nk.x(),nq*q_bk_nk.y(),nq*q_bk_nk.z());
     }
     ins.q_b_n=q_bk_nk;
     ins.C_b_n=Rotation::quaternion2rotationMatrix(q_bk_nk);
     return ins;
  }

















