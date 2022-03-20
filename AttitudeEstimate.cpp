/**
  ******************************************************************************
  * @file    AttitudeEstimate.CPP
  * @author  疯狂的兔子工作室
  * @version V1.0
  * @date    2022-3-8
  * @brief   姿态估计。
  ******************************************************************************
*/
#include "AttitudeEstimate.h"
#include "arm_math.h"
//#include "assert.h"

void AttitudeEstimate::Init_Matrix()//初始化矩阵尺寸
{
  //矩阵初始化（除Z V H HL HR1四个矩阵外，其他矩阵均初始化，H HL HR1在UpdateH中初始化，Z V在测量函数中初始化）
  arm_mat_init_f32(&A, 7, 7, A_data);//初始化A(7x7)矩阵（状态转移矩阵）
  arm_mat_init_f32(&UD_Q, 7, 1, UD_Q_data);
  arm_mat_init_f32(&UD_KAMA, 7, 7, UD_KAMA_data);
  arm_mat_init_f32(&UD_U, 7, 7, UD_U_data);
  arm_mat_init_f32(&X, 7, 1, X_data);//初始化状态向量(7x1)
  arm_mat_init_f32(&UD_D, 7, 1, UD_D_data);
  arm_mat_init_f32(&UD_R, 6, 1, UD_R_data);

  
  arm_mat_init_f32(&UD_Y, 7, 14, UD_Y_data);
  arm_mat_init_f32(&UD_Pnn, 7, 7, UD_Pnn_data);
  arm_mat_init_f32(&UD_D1, 14, 1, UD_D1_data);
  arm_mat_init_f32(&UD_V2N_1, 14, 1, UD_V2N_1_data);
  arm_mat_init_f32(&UD_V2N_2, 14, 1, UD_V2N_2_data);

  arm_mat_init_f32(&UD_ALPHA, 8, 1, UD_ALPHA_data);
  arm_mat_init_f32(&UD_ALPHA_L, 8, 1, UD_ALPHA_L_data);
  arm_mat_init_f32(&UD_F, 1, 7, UD_F_data);
  arm_mat_init_f32(&UD_V, 7, 1, UD_V_data);
  arm_mat_init_f32(&UD_FV, 7, 1, UD_FV_data);
  arm_mat_init_f32(&UD_B, 7, 1, UD_B_data);
  arm_mat_init_f32(&UD_M, 7, 1, UD_M_data);
  arm_mat_init_f32(&UD_Ki, 7, 1, UD_Ki_data);
  
  arm_mat_init_f32(&HR, 4, 3, HR_data);//初始化H矩阵辅助计算矩阵HR(4x3)
}
/**
  ******************************************************************************
  *   函 数 名: PreProcess(const arm_matrix_instance_f32* v_iPara,const float32_t* deltaT,const uint8_t* bAvalible)
  *   功能说明: 执行滤波前预处理，将输入参数做规范化处理，并放置到相应成员变量中
  *   形    参：v_iPara[0-2]：角速度测量值（顺序：X Y Z）
  *             v_iPara[3-5]：重力加速度测量值（顺序：X Y Z）
  *             v_iPara[6-8]：磁力计测量值（顺序：X Y Z）
  *             deltaT[0-8]：表示v_iPara[0-8]位置的参数与上一次滤波的时间差
  *             bAvalible[0-8]：表示v_iPara[0-8]位置的参数数据是否有效
  *   返 回 值: 无
  *   备 注：
  ********************************************************************************   
*/
void AttitudeEstimate::PreProcess(const arm_matrix_instance_f32* v_iPara,const float32_t* deltaT,const uint8_t* bAvalible)//
{
  //中心差分求时间步内的平均角速度
  arm_copy_f32(v_iPara->pData,fGyro,3);
  arm_add_f32(m_Gyro,fGyro,m_Gyro,3);//上一步角速度+本步角速度测量值
  arm_sub_f32(m_Gyro,X.pData+4,m_Gyro,3);//减去测量值偏差
  arm_scale_f32(m_Gyro,0.5f,m_Gyro,3);//除以2
  arm_copy_f32(deltaT,dt,9);
  //将加速度和地磁数据归一化
  if(bAvalible[3]&&bAvalible[4]&&bAvalible[5]){
    Normal(v_iPara->pData+4, fAcc, 3);
    m_fAcc_Available=1;
  }else{
    m_fAcc_Available=0;
  }
  if(bAvalible[6]&&bAvalible[7]&&bAvalible[8]){
    Normal(v_iPara->pData+7, fMag, 3);
    m_fMag_Available=1;
  }else{
    m_fMag_Available=0;
  }
}

void AttitudeEstimate::PostProcess()//执行滤波后处理
{
  //单位化q向量
  float32_t norm;
  arm_dot_prod_f32(X.pData,X.pData,4,&norm);
  arm_sqrt_f32(norm,&norm);
  norm=1.0f/norm;
  arm_scale_f32(X.pData,norm,X.pData,4);
  //根据角速度偏差更新角速度
  arm_sub_f32(fGyro,X.pData+4,m_Gyro,3);
  //更新四元数
  arm_copy_f32(X.pData,Quat.Q.pData,4);
}
/**
  ******************************************************************************
  *   函 数 名: Init_A()
  *   功能说明: 初始化A矩阵
  *   形    参：
  *   返 回 值: 无
  *   备 注：由于A矩阵为稀疏矩阵，本函数初始化A的单位矩阵和0矩阵分块部分。
  ********************************************************************************   
*/
void AttitudeEstimate::Init_A()//初始化A矩阵
{
  for(uint8_t i=0;i<49;i++){
    A_data[i] = 0.0f;
  }
  for(uint8_t i=0;i<7;i++){
    A_data[i*7+i] = 1.0f;
  }
}
/**
  ******************************************************************************
  *   函 数 名: Update_A()
  *   功能说明: 根据上一步的后验估计值初始化A矩阵
  *   形    参：无
  *   返 回 值: 无
  *   备 注：由于A矩阵为稀疏矩阵，函数只更新非0和非单位矩阵部分。
  ********************************************************************************   
*/
void AttitudeEstimate::Update_A()//计算A矩阵
{

  //计算A(7x7)矩阵
  float halfdt=dt/2.0f;

  //A(1:4,1:4)=I+OMIGA;
  //           |0   -Ox  -Oy  -Oz|
  //OMIGA=dt/2*|Ox   0    Oz  -Oy|
  //           |Oy  -Oz   0    Ox|
  //           |Oz   Oy	 -Ox   0 |
  //I为单位矩阵
  float Ox=m_Gyro[0]*halfdt;
  float Oy=m_Gyro[1]*halfdt;
  float Oz=m_Gyro[2]*halfdt;
  A.pData[0]=1.0f;  A.pData[1]=-Ox;     A.pData[2]=-Oy;     A.pData[3]=-Oz;
  A.pData[7]=Ox;    A.pData[8]=1.0f;    A.pData[9]=Oz;      A.pData[10]=-Oy;
  A.pData[14]=Oy;   A.pData[15]=-Oz;    A.pData[16]=1.0f;   A.pData[17]=Ox;
  A.pData[21]=Oz;   A.pData[22]=Oy;     A.pData[23]=-Ox;    A.pData[24]=1.0f;

  float q0hafT = X.pData[0] * halfdt;
  float q1hafT = X.pData[1] * halfdt;
  float q2hafT = X.pData[2] * halfdt;
  float q3hafT = X.pData[3] * halfdt;
  
  //                | q1   q2   q3|
  //A(1:4,5:7)=dt/2*|-q0   q3  -q2|
  //                |-q3  -q0   q1|
  //                | q2  -q1  -q0|
  A.pData[4] =   q1hafT; A.pData[5] =   q2hafT; A.pData[6] =   q3hafT;
  A.pData[11] = -q0hafT; A.pData[12] =  q3hafT; A.pData[13] = -q2hafT;
  A.pData[18] = -q3hafT; A.pData[19] = -q0hafT; A.pData[20] =  q1hafT;
  A.pData[25] =  q2hafT; A.pData[26] = -q1hafT; A.pData[27] = -q0hafT;
  //A(5:7,1:7)=0(1:3,1:4)|I(1:3,1:3)，初始化时进行了赋值

  //计算H的后三列(H后三列使用状态转移前的参数计算)
  //                | q1   q2   q3|
  //HR(4x3)= dt/2 * |-q0   q3  -q2|
  //                |-q3  -q0   q1|
  //                | q2  -q1  -q0|
  HR.pData[0] =  q1hafT; HR.pData[1] =  q2hafT; HR.pData[2] =  q3hafT;
  HR.pData[3] = -q0hafT; HR.pData[4] =  q3hafT; HR.pData[5] = -q2hafT;
  HR.pData[6] = -q3hafT; HR.pData[7] = -q0hafT; HR.pData[8] =  q1hafT;
  HR.pData[9] =  q2hafT; HR.pData[10] = -q1hafT; HR.pData[11] = -q0hafT;
}

/**
  ******************************************************************************
  *   函 数 名: Update_H()
  *   功能说明: 根据上一步的后验估计值初始化H矩阵
  *   形    参：
  *   返 回 值: 无
  *   备 注：
  ********************************************************************************   
*/

void AttitudeEstimate::Update_H()
{
  float hx,hy,by,bz;
  
  float q0=X.pData[0];
  float q1=X.pData[1];
  float q2=X.pData[2];
  float q3=X.pData[3];
  if(m_fMag_Available){
    arm_mat_init_f32(&H, 6, 7, H_data);//初始化H矩阵
    arm_mat_init_f32(&HL, 6, 4, HL_data);//初始化H矩阵辅助计算矩阵HL
    arm_mat_init_f32(&HR1, 6, 3, HR1_data);//初始化H矩阵辅助计算矩阵HR1


    float q0q1 = q0*q1;
    float q0q2 = q0*q2;
    float q0q3 = q0*q3;
    float q1q1 = q1*q1;
    float q1q2 = q1*q2;
    float q1q3 = q1*q3;
    float q2q2 = q2*q2;
    float q2q3 = q2*q3;
    float q3q3 = q3*q3;
    
    //将体坐标系下测量的地磁场强度数据转换到参考坐标系下
    hx=2.0f*(fMag[0]*(0.5f-q2q2-q3q3)+fMag[1]*(q1q2-q0q3)+fMag[2]*(q1q3+q0q2));
    hy=2.0f*(fMag[0]*(q1q2+q0q3)+fMag[1]*(0.5f-q1q1-q3q3)+fMag[2]*(q2q3-q0q1));
    bz=2.0f*(fMag[0]*(q1q3-q0q2)+fMag[1]*(q2q3+q0q1)+fMag[2]*(0.5f-q1q1-q2q2));
    //对参考坐标系下的磁场强度进行修正
    arm_sqrt_f32(hx*hx+hy*hy,&by);
    float q0by=q0*by; float q1by=q1*by; float q2by=q2*by; float q3by=q3*by;
    float q0bz=q0*bz; float q1bz=q1*bz; float q2bz=q2*bz; float q3bz=q3*bz;
    
    HL.pData[0] = -q2;  HL.pData[1] =  q3;  HL.pData[2] = -q0;  HL.pData[3] = q1;
    HL.pData[4] =  q1;  HL.pData[5] =  q0;  HL.pData[6] =  q3;  HL.pData[7] = q2;
    HL.pData[8] =  q0;  HL.pData[9] = -q1;  HL.pData[10] = -q2;  HL.pData[11] = q3;
    
    HL.pData[12]= q3by-q2bz;  HL.pData[13]= q2by+q3bz;  HL.pData[14]=q1by-q0bz;  HL.pData[15]= q0by+q1bz;
    HL.pData[16]= q0by+q1bz;  HL.pData[17]=-q1by+q0bz;  HL.pData[18]=q2by+q3bz;  HL.pData[19]=-q3by+q2bz;
    HL.pData[20]=-q1by+q0bz;  HL.pData[21]=-q0by-q1bz;  HL.pData[22]=q3by-q2bz;  HL.pData[23]= q2by+q3bz;
    //HL *= 2.0f;
    for(uint8_t i=0;i<HL.numRows*HL.numCols;i++){
      HL.pData[i]*=2.0f;
    }
    //HR1 = HL * HR;
    arm_mat_mult_f32(&HL, &HR, &HR1);
    //合并进H矩阵
    for (uint8_t i = 0; i < H.numRows; i++) {
      for (uint8_t j = 0; j < HL.numCols; j++) {
        H.pData[i*H.numCols+j] = HL.pData[i*HL.numCols+j];
      }
      for (uint8_t j = 0; j < HR1.numCols; j++) {
        H.pData[i*H.numCols+(j+HL.numCols)] = HR1.pData[i*HR1.numCols+j];
      }
    }
  }else{
    arm_mat_init_f32(&H, 3, 7, H_data);//初始化H矩阵
    arm_mat_init_f32(&HL, 3, 4, HL_data);//初始化H矩阵辅助计算矩阵HL
    arm_mat_init_f32(&HR1, 3, 3, HR1_data);//初始化H矩阵辅助计算矩阵HR1

    HL.pData[0] = -q2;  HL.pData[1] =  q3;  HL.pData[2] = -q0;  HL.pData[3] = q1;
    HL.pData[4] =  q1;  HL.pData[5] =  q0;  HL.pData[6] =  q3;  HL.pData[7] = q2;
    HL.pData[8] =  q0;  HL.pData[9] = -q1;  HL.pData[10] = -q2;  HL.pData[11] = q3;
    //HL *= 2.0f;
    for(uint8_t i=0;i<HL.numRows*HL.numCols;i++){
      HL.pData[i]*=2.0f;
    }
    //HR1 = HL * HR;
    arm_mat_mult_f32(&HL, &HR, &HR1);
    //合并进H矩阵
    for (uint8_t i = 0; i < H.numRows; i++) {
      for (uint8_t j = 0; j < HL.numCols; j++) {
        H.pData[i*H.numCols+j] = HL.pData[i*HL.numCols+j];
      }
      for (uint8_t j = 0; j < HR1.numCols; j++) {
        H.pData[i*H.numCols+(j+HL.numCols)] = HR1.pData[i*HR1.numCols+j];
      }
    }
  }
}
/**
  ******************************************************************************
  *   函 数 名: Init_P()
  *   功能说明: 初始化协方差矩阵P=UDUT
  *   形    参：
  *   返 回 值: 无
  *   备 注：由于卡尔曼滤波算法是迭代算法，P矩阵最终收敛在真实值，因此P矩阵的初始值对结果无影响。
  *          但P矩阵不能初始化为奇异矩阵，即|P|!=0
  ********************************************************************************   
*/
void AttitudeEstimate::Init_P()
{
  memset(UD_U.pData,0,sizeof(float32_t)*UD_U.numCols*UD_U.numRows);//清零
  memset(UD_KAMA.pData,0,sizeof(float32_t)*UD_KAMA.numCols*UD_KAMA.numRows);//清零
  for(uint16_t i=0;i<UD_U.numRows;i++){
    UD_U.pData[i*UD_U.numCols+i] = 1.0f;
    UD_D.pData[i]=1.0f;
    UD_KAMA.pData[i*UD_KAMA.numCols+i] = 1.0f;
  }
}

/**
  ******************************************************************************
  *   函 数 名: Init_Q()
  *   功能说明: 初始化系统噪声协方差矩阵
  *   形    参：无
  *   返 回 值: 无
  *   备 注：系统噪声协方差矩阵对整体计算影响不大，但Q越准确，滤波估计就越准确
  ********************************************************************************   
*/
void AttitudeEstimate::Init_Q()
{
  for(uint8_t i=0;i<4;i++){
    UD_Q.pData[i]=1.0e-6;
  }
  for(uint8_t i=4;i<7;i++){
    UD_Q.pData[i]=1.0e-10;
  }
}
/**
  ******************************************************************************
  *   函 数 名: Init_R(const arm_matrix_instance_f32* v_iPara)
  *   功能说明: 初始化测量传感器协方差矩阵
  *   形    参：无
  *   返 回 值: 无
  *   备 注：测量噪声协方差矩阵对整体计算影响不大，但R越准确，滤波估计就越准确
  *          可以通过分析传感器的噪声获得测量噪声协方差矩阵
  ********************************************************************************   
*/
void AttitudeEstimate::Init_R()//初始化R矩阵
{
  UD_R.pData[0]=1.5f;UD_R.pData[1]=1.5f;UD_R.pData[2]=1.5f;//加速度计方差（假设各方向相互独立）
  UD_R.pData[3]=75.0f;UD_R.pData[4]=75.0f;UD_R.pData[5]=75.0f;//地磁传感器方差（假设各方向相互独立）
}
/**
  ******************************************************************************
  *   函 数 名: Init_X(const arm_matrix_instance_f32* v_iPara)
  *   功能说明: 初始化状态向量
  *   形    参：v_iPara 四元数向量   v_iPara.pData[i]=qi  i=0,1,2,3,4
  *   返 回 值: 无
  *   备 注：
  ********************************************************************************   
*/
void AttitudeEstimate::Init_X(const arm_matrix_instance_f32* v_iPara)//初始化状态矩阵
{
  arm_copy_f32(v_iPara->pData,X.pData,4);
  arm_copy_f32(v_iPara->pData,Quat.Q.pData,4);
  X.pData[4]=1e-9;
  X.pData[5]=1e-9;
  X.pData[6]=1e-9;
}
/**
  ******************************************************************************
  *   函 数 名: StateEquation()
  *   功能说明: 状态转移方程，利用上一步估计的四元数和角速度值，
  *             利用一阶龙格库塔法积分获得本次的四元数。
  *   形    参：无
  *   返 回 值: 无
  *   备 注：
  ********************************************************************************
*/
void AttitudeEstimate::StateEquation()
{
  float norm;
  float halfT = dt / 2.0f;
  float q0 = X.pData[0];
  float q1 = X.pData[1];
  float q2 = X.pData[2];
  float q3 = X.pData[3];
  //利用一阶龙格库塔法进行积分（此处采用矩形积分）
  //，获得姿态四元数，其中halfT为积分步长（即采样周期）的一半。					   
  q0 = q0 + (-q1 * m_Gyro[0] - q2 * m_Gyro[1] - q3 * m_Gyro[2])*halfT;
  q1 = q1 + (q0*m_Gyro[0] + q2 * m_Gyro[2] - q3 * m_Gyro[1])*halfT;
  q2 = q2 + (q0*m_Gyro[1] - q1 * m_Gyro[2] + q3 * m_Gyro[0])*halfT;
  q3 = q3 + (q0*m_Gyro[2] + q1 * m_Gyro[1] - q2 * m_Gyro[0])*halfT;
  arm_sqrt_f32(q0*q0 + q1 * q1 + q2 * q2 + q3 * q3,&norm);
  q0 = q0 / norm;
  q1 = q1 / norm;
  q2 = q2 / norm;
  q3 = q3 / norm;
  X.pData[0] = q0;
  X.pData[1] = q1;
  X.pData[2] = q2;
  X.pData[3] = q3;
}
/**
  ******************************************************************************
  *   函 数 名: MeasuringEquation()
  *   功能说明: 测量预测方程
  *   形    参：
  *   返 回 值: 无
  *   备 注：
  ********************************************************************************
*/
void AttitudeEstimate::MeasuringEquation()
{
  float hx, hy, by, bz;
  float q0 = X.pData[0];
  float q1 = X.pData[1];
  float q2 = X.pData[2];
  float q3 = X.pData[3];
  float q0q0 = q0 * q0;
  float q0q1 = q0 * q1;
  float q0q2 = q0 * q2;
  float q1q1 = q1 * q1;
  float q1q3 = q1 * q3;
  float q2q2 = q2 * q2;
  float q2q3 = q2 * q3;
  float q3q3 = q3 * q3;
  if (m_fMag_Available) {
    //形成测量向量实测值Z
    arm_mat_init_f32(&Z,6,1,Z_data);
    for (uint8_t i = 0; i < 3; i++) {
      Z.pData[i] = fAcc[i];
      Z.pData[i+3]=fMag[i];
    }
    //形成测量向量预测值V
    arm_mat_init_f32(&V,6,1,V_data);
    V.pData[0] = 2 * (q1q3 - q0q2);
    V.pData[1] = 2 * (q0q1 + q2q3);
    V.pData[2] = q0q0 - q1q1 - q2q2 + q3q3;
    float q0q3 = q0 * q3;
    float q1q2 = q1 * q2;
    //将体坐标系下测量的地磁场强度数据转换到参考坐标系下
    hx = 2.0f*(fMag[0] * (0.5f - q2q2 - q3q3) + fMag[1] * (q1q2 - q0q3) + fMag[2] * (q1q3 + q0q2));
    hy = 2.0f*(fMag[0] * (q1q2 + q0q3) + fMag[1] * (0.5f - q1q1 - q3q3) + fMag[2] * (q2q3 - q0q1));
    bz = 2.0f*(fMag[0] * (q1q3 - q0q2) + fMag[1] * (q2q3 + q0q1) + fMag[2] * (0.5f - q1q1 - q2q2));
    //对参考坐标系下的磁场强度进行修正
    arm_sqrt_f32(hx*hx + hy * hy,&by);
    //将修正后的磁场强度转换到机体坐标系下
    V.pData[3] = 2.0f*by*(q1q2 + q0q3) + 2.0f*bz*(q1q3 - q0q2);
    V.pData[4] = 2.0f*by*(0.5f - q1q1 - q3q3) + 2.0f*bz*(q0q1 + q2q3);
    V.pData[5] = 2.0f*by*(q2q3 - q0q1) + 2.0f*bz*(0.5f - q1q1 - q2q2);
  }else{
    //形成测量向量实测值Z
    arm_mat_init_f32(&Z,3,1,Z_data);
    for (uint8_t i = 0; i < 3; i++) {
      Z.pData[i] = fAcc[i];
    }
    //形成测量向量预测值V
    arm_mat_init_f32(&V,3,1,V_data);
    V.pData[0] = 2 * (q1q3 - q0q2);
    V.pData[1] = 2 * (q0q1 + q2q3);
    V.pData[2] = q0q0 - q1q1 - q2q2 + q3q3;
  }
}
/**
  ******************************************************************************
  *   函 数 名: AttitudeEstimate()
  *   功能说明: 构造函数
  *   形    参：无
  *   返 回 值: 无
  *   备 注：
  ********************************************************************************   
*/
AttitudeEstimate::AttitudeEstimate()
{
  for (uint8_t i = 0; i < 3; i++) {
    m_Gyro[i] = 0.0f;
  }
}
/**
  ******************************************************************************
  *   函 数 名: ~AttitudeEstimate()
  *   功能说明: 构造函数
  *   形    参：无
  *   返 回 值: 无
  *   备 注：
  ********************************************************************************   
*/
AttitudeEstimate::~AttitudeEstimate()
{
}

/**
  ******************************************************************************
  *   函 数 名: Normal(float* fromV,float* toV, uint8_t len)
  *   功能说明: 向量归一化
  *   形    参：float* fromV被归一化的向量，float* toV归一化后的向量，uint8_t len向量长度
  *   返 回 值: 无
  *   备 注：
  ********************************************************************************
*/

void AttitudeEstimate::Normal(float* fromV,float* toV, uint8_t len)
{
  float32_t norm;
  arm_dot_prod_f32(fromV,fromV,len,&norm);
  arm_sqrt_f32(norm,&norm);
  if(norm>0.0f){//如果向量的长度大于0，则进行单位化，否则不对向量做处理
    norm=1.0f/norm;
    arm_scale_f32(fromV,norm,toV,len);
  }
}
