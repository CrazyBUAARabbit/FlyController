/**
  ******************************************************************************
  * @file    AttitudeEstimate.h
  * @author  疯狂的兔子工作室
  * @version V1.0
  * @date    2022-3-8
  * @brief   姿态估计。
  ******************************************************************************
  */
#ifndef _AttitudeEstimate_H_
#define _AttitudeEstimate_H_
#include "stm32f3xx_hal.h"
#include "CKalmanFilter.h"
#include "Quaternion.h"
#include "arm_math.h"
class AttitudeEstimate:public CKalmanFilter
{
private:
  void Normal(float* fromV, float* toV, uint8_t len);
protected:
  //主要数据区1(共计:228*4Byte)
  float32_t A_data[49],UD_Q_data[7],UD_KAMA_data[49],UD_U_data[49];
  float32_t H_data[42],X_data[7],UD_D_data[7],Z_data[6],V_data[6],UD_R_data[6];
  
  //卡尔曼滤波计算过程辅助数据区(共计:233*4Byte)
  float32_t UD_Y_data[98],UD_Pnn_data[49],UD_D1_data[14],UD_V2N_1_data[14],UD_V2N_2_data[14];
  float32_t UD_ALPHA_data[8],UD_ALPHA_L_data[8],UD_F_data[7],UD_V_data[7],UD_FV_data[7];
  float32_t UD_B_data[7],UD_M_data[7],UD_Ki_data[7];
  
  //H矩阵计算辅助数据区(共计:54*4Byte)
  arm_matrix_instance_f32 HR;
  float32_t HR_data[12];
  arm_matrix_instance_f32 HL,HR1;
  float32_t HL_data[24],HR1_data[18];
public:
  AttitudeEstimate();
  ~AttitudeEstimate();
  Quaternion Quat;//四元数
  float32_t m_Gyro[3];//当前角速度（滤波后）

protected:
  float32_t fGyro[3];//角速度传感器测量值（滤波前）
  float32_t dt[9];//距离上次计算的时间差dt
  float32_t fAcc[3];//加速度传感器测量值（滤波前）
  float32_t fMag[3];//地磁场传感器测量值（滤波前）
  uint8_t m_fMag_Available;//地磁传感器测量值是否有效标志位  0-无效   非0-有效
  uint8_t m_fAcc_Available;//加速度传感器测量值是否有效标志位  0-无效   非0-有效
protected:
  virtual void Init_Matrix();//初始化矩阵尺寸
  virtual void PreProcess(const arm_matrix_instance_f32* v_iPara,const float32_t* deltaT,const uint8_t* bAvalible);//执行滤波前预处理
  virtual void Init_A();//初始化状态转移A矩阵
  virtual void Update_A();//更新状态转移A矩阵
  //virtual void Init_H();//初始化测量转换H矩阵
  virtual void Update_H();//更新测量转换H矩阵
  virtual void Init_P();//初始化状态向量方差P矩阵
  virtual void Init_Q();//初始化过程噪声方差Q矩阵
  //virtual void Update_Q();//更新过程噪声方差Q矩阵
  virtual void Init_R();//初始化测量噪声方差R矩阵（假设测量向量相互独立，即R为对角阵）
  virtual void Update_R();//更新测量噪声方差R矩阵（假设测量向量相互独立，即R为对角阵）
  virtual void Init_X(const arm_matrix_instance_f32* v_iPara);//初始化状态向量X
  virtual void StateEquation();//状态转移方程
  virtual void MeasuringEquation();//测量预测方程
  virtual void PostProcess();//执行滤波后处理
};
#endif