/**
  ******************************************************************************
  * @file    CKalmanFilter.h
  * @author  疯狂的兔子工作室
  * @version V1.1
  * @date    2022-1-16
  * @brief   扩展卡尔曼滤波器.
  ******************************************************************************
  */
#ifndef _CKalmanFilter_H_
#define _CKalmanFilter_H_
#include "stm32f3xx_hal.h"
#include "arm_math.h"
class CKalmanFilter
{
protected:
  //主要数据区
  arm_matrix_instance_f32 A;//状态转移A矩阵(nxn)
  arm_matrix_instance_f32 UD_Q;//过程噪声方差Q矩阵，Q矩阵为对角矩阵，当X变量相关时，变换KAMA矩阵，使Q转化为对角阵，(nx1)
  arm_matrix_instance_f32 UD_KAMA;//KAMA矩阵(nxn)
  arm_matrix_instance_f32 UD_U;//协方差矩阵P的UD分解值(nxn)
  arm_matrix_instance_f32 H;//测量转换H矩阵(mxn)
  arm_matrix_instance_f32 X;//状态向量(nx1)
  arm_matrix_instance_f32 UD_D;//协方差矩阵P的UD分解值(nx1)
  arm_matrix_instance_f32 Z;//传感器实际测量值（向量）(mx1)
  arm_matrix_instance_f32 V;//传感器预测测量值（向量）(mx1)
  arm_matrix_instance_f32 UD_R;//测量噪声方差R矩阵（假设测量向量相互独立，即R为对角阵），R为对角线元素(mx1)
  
  //卡尔曼滤波计算过程辅助数据区
  arm_matrix_instance_f32 UD_Y;//(nx2n)
  arm_matrix_instance_f32 UD_Pnn;//(nxn)
  arm_matrix_instance_f32 UD_D1;//(2nx1)
  arm_matrix_instance_f32 UD_V2N_1;//(2nx1)
  arm_matrix_instance_f32 UD_V2N_2;//(2nx1)
  arm_matrix_instance_f32 UD_ALPHA;//((n+1)x1)
  arm_matrix_instance_f32 UD_ALPHA_L;//((n+1)x1)
  arm_matrix_instance_f32 UD_F;//(nx1)
  arm_matrix_instance_f32 UD_V;//(nx1)
  arm_matrix_instance_f32 UD_FV;//(nx1)
  arm_matrix_instance_f32 UD_B;//(nx1)
  arm_matrix_instance_f32 UD_M;//(nx1)
  arm_matrix_instance_f32 UD_Ki;//(nx1)

  virtual void Init_Matrix();//初始化矩阵尺寸
  virtual void PreProcess(const arm_matrix_instance_f32* v_iPara,const float32_t* deltaT,const uint8_t* bAvalible);//执行滤波前预处理
  virtual void PostProcess();//执行滤波后处理
  virtual void Init_A();//初始化状态转移A矩阵
  virtual void Update_A();//更新状态转移A矩阵
  virtual void Init_H();//初始化测量转换H矩阵
  virtual void Update_H();//更新测量转换H矩阵
  virtual void Init_P();//初始化状态向量方差P矩阵
  virtual void Init_Q();//初始化过程噪声方差Q矩阵
  virtual void Update_Q();//更新过程噪声方差Q矩阵
  virtual void Init_R();//初始化测量噪声方差R矩阵（假设测量向量相互独立，即R为对角阵）
  virtual void Update_R();//更新测量噪声方差R矩阵（假设测量向量相互独立，即R为对角阵）
  virtual void Init_X(const arm_matrix_instance_f32* v_iPara);//初始化状态向量X
  virtual void StateEquation();//状态转移方程
  virtual void MeasuringEquation();//测量预测方程
public:
  CKalmanFilter();
  ~CKalmanFilter();
  void Init(const arm_matrix_instance_f32* v_iPara);//对卡尔曼滤波器进行初始化
  void FilterExe(const arm_matrix_instance_f32* v_iPara,const float32_t* deltaT,const uint8_t* bAvalible);//执行卡尔曼滤波
};
#endif