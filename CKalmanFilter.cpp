/**
  ******************************************************************************
  * @file    CKalmanFilter.CPP
  * @author  疯狂的兔子工作室
  * @version V1.1
  * @date    2022-3-6
  * @brief   扩展卡尔曼滤波。
  ******************************************************************************
*/
/*
版本说明：
V1.1版与V1.0版区别：
1.将原滤波函数更改为基于序贯UD分解卡尔曼滤波算法，使滤波算法成为通用算法；
2.将卡尔曼滤波算法改为基类，需要使用卡尔曼滤波时，只需要继承该基类即可；
3.为达到效果，将存储空间放在子类中定义将初始化函数、更新函数设置为虚函数，将各矩阵的内存空间声明放入子类。
*/
#include "CKalmanFilter.h"
#include "arm_math.h"
#include "assert.h"
#define NDEBUG

CKalmanFilter::CKalmanFilter()
{
}
CKalmanFilter::~CKalmanFilter()
{
}
/**
  ******************************************************************************
  *   函 数 名: FilterExe
  *   功能说明: 基于UD分解的序贯扩展的卡尔曼滤波器，较标志的扩展卡尔曼滤波器计算量较小(序贯滤波优势)，且X的自协方差矩阵P以U、D存储，
  *             解决因算法舍入误差导致P非对称或非正定问题
  *   形    参：滤波器输入值向量:
                v_iPara向量的具体定义由子类成员函数PreProcess的形参明确
                deltaT 向量表示v_iPara向量中的物理量的时间步长（即与上一次滤波时的时间差）
                bAvalible表示本次v_iPara向量相应分量的物理量的有效性，1-表示有效，0-表示无效
  *   返 回 值: 无
  *   备 注：状态向量为四元数、陀螺仪偏差
  *          计算步骤为：
  *          0.根据状态方程f(X)对X上一步的X无偏估计值求偏微分更新A矩阵
  *          1.计算预测值与真实值之间的协方差矩阵UD分解值（先验估计）：P=A*P*AT+Q=A*U*D*UT*AT+Q=Y*D1*YT
  *            Y=|Ak*Uk-1  I|  D1=|Dk-1  0|
  *                               |0     Q|
  *            根据Gram-Schmidt公式计算出Dk/k-1，Uk/k-1
  *            for j=n,n-1,...,1
  *             Dk/k-1(j,j)=Sigma(D1(s,s)*Y(j,s)^2  s=1,2,...,2n)   其中n为Ak的维度
  *               for i=1,2,...,j-1
  *                 Uk/k-1(i,j)=Sigma(D1(s,s)*Y(i,s)*Y(j,s)  s=1,2,...,2n)/(Dk/k-1(j,j))
  *                 Y(i,:)=Y(i,:)-Uk/k-1(i,j)*Y(j,:)
  *               end i
  *             end j
  *          2.计算状态向量（先验估计）：X=f(X)
  *          3.根据测量方程g(X)在X先验估计值处的偏微分更新H矩阵
  *          4.对测量矩阵H进行序贯处理，即取出H的行向量H(i)，进行m次迭代对X和UD进行测量修正(后验估计)
  *            注：序贯处理的前提条件为R矩阵为对角阵，即测量向量个变量相互独立。
  *          5.计算协方差矩阵Uk,Dk（后验估计）:
                F=UT*H(i)T=(H(i)*U)T，其中H(i)为H的第i行；
  *             V=Dk/k-1*F  得出V(j)=Dk/k-1(j,j)*F(j)   j=1,2,...m
  *             ALPHA(0)=R(i,i)
  *             for l=1,2,3,...n
  *               ALPHA(l)=ALPHA(l-1)+F(l)*V(l)
  *               B(l)=V(l)
  *               M(l)=-F(l)/ALPHA(l-1)
  *               Dk(l,l)=Dk/k-1(l,l)*ALPHA(l-1)/ALPHA(l)
  *               for o=1,2,...,l-1
  *                  Uk(o,l)=Uk/k-1(o,l)+B(o)*M(l)
  *                  B(o)+=Uk(o,l)*V(l)
  *               end o
  *             end l
  *          6.计算卡尔曼增益：K=B/ALPHA(n)
  *          7.计算状态向量（后验估计）：X=X+K*(Z(i)-g(X,i))
  ********************************************************************************   
*/
void CKalmanFilter::FilterExe(const arm_matrix_instance_f32* v_iPara,const float32_t* deltaT,const uint8_t* bAvalible)
{
  
  //判定主数据区尺寸约束
  assert(A.numCols==A.numRows);//判定A为方阵，(nxn)
  assert(UD_U.numCols==A.numCols&&UD_U.numRows==A.numRows);//判定A、U尺寸相等，(nxn)
  assert(A.numRows==UD_KAMA.numRows&&A.numCols==UD_KAMA.numCols);//判定A、KAMA尺寸相等，(nxn)
  assert(UD_Q.numCols==1&&UD_Q.numRows==A.numRows);//判定Q为列向量，且行数与A、P相等(nx1)
  assert(X.numRows==A.numRows&&X.numCols==1);//判定X为列向量，且行数与A、P相等(nx1)
  assert(UD_D.numCols==1&&UD_D.numRows==A.numRows);//D为列向量，长度与A的边长相等(nx1)
  //判定计算辅助区尺寸约束
  assert(UD_Y.numCols==A.numCols*2&&UD_Y.numRows==A.numRows);//(nx2n)
  assert(UD_Pnn.numCols==A.numCols&&UD_Pnn.numRows==A.numRows);//(nxn)
  assert(UD_D1.numRows==A.numRows*2&&UD_D1.numCols==1);//(2nx1)
  assert(UD_V2N_1.numRows==A.numRows*2&&UD_V2N_1.numCols==1);//(2nx1)
  assert(UD_V2N_2.numRows==A.numRows*2&&UD_V2N_2.numCols==1);//(2nx1)
  assert(UD_ALPHA.numRows==A.numRows+1&&UD_ALPHA.numCols==1);//((n+1)x1)
  assert(UD_ALPHA_L.numRows==A.numRows+1&&UD_ALPHA_L.numCols==1);//((n+1)x1)
  assert(UD_F.numRows==1&&UD_F.numCols==A.numRows);//(1xn)
  assert(UD_V.numRows==A.numRows&&UD_V.numCols==1);//(nx1)
  assert(UD_FV.numRows==A.numRows&&UD_FV.numCols==1);//(nx1)
  assert(UD_B.numRows==A.numRows&&UD_B.numCols==1);//(nx1)
  assert(UD_M.numRows==A.numRows&&UD_M.numCols==1);//(nx1)
  assert(UD_Ki.numRows==A.numRows&&UD_Ki.numCols==1);//(nx1)
  
  PreProcess(v_iPara,deltaT,bAvalible);//对输入参数进行预处理
  
  /**** 1.计算协方差矩阵P（用U、D替代）和状态转移（先验估计） ****/
  /**** 1.1 U、D先验估计 ****/
  //更新A矩阵
  Update_A();
  //更新系统噪声协方差矩阵
  Update_Q();
  //构建Y=|Ak*Uk-1  KAMA|，D1=|D(对角元素) Q(对角元素)|
  arm_mat_mult_f32(&A, &UD_U, &UD_Pnn);
  for(uint16_t i=0;i<A.numRows;i++){
    arm_copy_f32(&(UD_Pnn.pData[i*A.numRows]),&(UD_Y.pData[i*UD_Y.numCols]),A.numRows);
    arm_copy_f32(&(UD_KAMA.pData[i*A.numRows]),&(UD_Y.pData[i*UD_Y.numCols+A.numRows]),A.numRows);
    UD_D1.pData[i]=UD_D.pData[i];//D使用列向量存储对角线
    UD_D1.pData[i+A.numRows]=UD_Q.pData[i];//Q使用列向量存储对角线
  }
  //根据Gram-Schmidt公式计算出Dk/k-1，Uk/k-1
  for(int16_t j=A.numRows-1;j>=0;j--){//for j=n,n-1,...,1
    //c(m)=D1(m,m)*Y(j,m)   m=1,2,...,2n   其中n为Ak的维度
    arm_mult_f32(UD_D1.pData,&(UD_Y.pData[j*UD_Y.numCols]),UD_V2N_1.pData,UD_Y.numCols);
    //Dk/k-1(j,j)=Sigma(c(s)*Y(j,s)  s=1,2,...,2n)
    arm_dot_prod_f32(UD_V2N_1.pData,&(UD_Y.pData[j*UD_Y.numCols]),UD_Y.numCols,&(UD_D.pData[j]));
    //c(m)=c(m)/Dk/k-1(j,j)   m=1,2,...,2n
    arm_scale_f32(UD_V2N_1.pData,1.0f/UD_D.pData[j],UD_V2N_1.pData,UD_Y.numCols);
                  
    for(uint16_t i=0;i<j-1;i++){//for i=1,2,...,j-1
      //Uk/k-1(i,j)=Sigma(c(s)*Y(i,s)  s=1,2,...,2n)
      arm_dot_prod_f32(UD_V2N_1.pData,&(UD_Y.pData[i*UD_Y.numCols]),UD_Y.numCols,&(UD_U.pData[i*A.numRows+j]));
      
      //Y(i,:)=Y(i,:)-Uk/k-1(i,j)*Y(j,:)
      arm_scale_f32(&(UD_Y.pData[j*UD_Y.numCols]),UD_U.pData[i*A.numRows+j],UD_V2N_2.pData,UD_Y.numCols);
      arm_sub_f32(&(UD_Y.pData[i*UD_Y.numCols]),UD_V2N_2.pData,&(UD_Y.pData[i*UD_Y.numCols]),UD_Y.numCols);
    }
  }
  
  /****~1.1 U、D先验估计结束****/
  /****1.2 状态转移 ****/
  StateEquation();
  /****~1.2 状态转移结束  **/
  /**** ~1.先验估计计算结束 ****/
  
  //更新H矩阵
  Update_H();
  //计算测量值Z、测量预测值向量Z
  MeasuringEquation();
  
  assert(Z.numRows==V.numRows&&V.numCols==Z.numCols);//判定V与Z尺寸相等(mx1)

  
  /****2.使用序贯法进行后验估计****/
  arm_matrix_instance_f32 UD_Hi;//辅助计算数据
  for(uint16_t i=0;i<H.numRows;i++){//使用序贯方法逐行处理H:H(i)
    /*******2.1 数据准备**********/
    arm_mat_init_f32(&UD_Hi, 1, 7, &(H.pData[i*A.numRows]));//取出H的第i行，作为向量UD_Hi
    //F=(H(i)*U)T
    arm_mat_mult_f32(&UD_Hi, &UD_U, &UD_F);
    //V(j)=Dk/k-1(j,j)*F(j)   j=1,2,...m
    arm_mult_f32(UD_D.pData,UD_F.pData,UD_V.pData,A.numRows);
    //ALPHA(0)=R(i,i)
    UD_ALPHA.pData[0]=UD_R.pData[i];
    UD_ALPHA_L.pData[0]=1.0f/UD_ALPHA.pData[0];
    //ALPHA(l)=ALPHA(l-1)+F(l)*V(l)   l=1,2,3,...n
    arm_mult_f32(UD_F.pData,UD_V.pData,UD_FV.pData,A.numRows);
    for(uint16_t l=0;l<A.numRows;l++){
      UD_ALPHA.pData[l+1]=UD_ALPHA.pData[l]+UD_FV.pData[l];
      UD_ALPHA_L.pData[l+1]=1.0f/UD_ALPHA.pData[l+1];//UD_ALPHA_L中保存UD_ALPHA的倒数
    }
    //B(l)=V(l)
    arm_copy_f32(UD_V.pData,UD_B.pData,A.numRows);
    //M(l)=-F(l)/ALPHA(l-1)
    arm_mult_f32(UD_F.pData,UD_ALPHA_L.pData,UD_M.pData,A.numRows);
    arm_negate_f32(UD_M.pData,UD_M.pData,A.numRows);
    /*******~2.1 数据准备完毕**********/
    
    /*******2.2 U、D后验估计**********/
    //Dk(l,l)=Dk/k-1(l,l)*ALPHA(l-1)/ALPHA(l)
    arm_mult_f32(UD_D.pData,UD_ALPHA.pData,UD_D.pData,A.numRows);
    arm_mult_f32(UD_D.pData,UD_ALPHA_L.pData+1,UD_D.pData,A.numRows);
    
    for(uint16_t l=0;l<A.numRows;l++){////for l=1,2,3,...n
      for(uint16_t o=0;o<l-1;o++){//for o=1,2,...,l-1
        //Uk(o,l)=Uk/k-1(o,l)+B(o)*M(l)
        UD_U.pData[o*A.numRows+l]+=UD_B.pData[o]*UD_M.pData[l];
        //B(o)+=Uk(o,l)*V(l)
        UD_B.pData[o]+=UD_U.pData[o*A.numRows+l]*UD_V.pData[l];
      }
    }
    /*******~2.2 U、D后验估计完毕**********/
    
    /********2.3计算卡尔曼增益*******/
    //K=B/ALPHA(n)
    arm_scale_f32(UD_B.pData,UD_ALPHA_L.pData[A.numRows],UD_Ki.pData,A.numRows);
    /********~2.3计算卡尔曼增益结束*******/
    
    /****2.4 X=X+K*(Z-H*X)状态向量后验估计****/  
    //Z=Z-g(X)
    Z.pData[i]-= V.pData[i];
    //K = K * Z;
    arm_scale_f32(UD_Ki.pData,Z.pData[i],UD_Ki.pData,A.numRows);
    //X += K;
    arm_add_f32(X.pData,UD_Ki.pData,X.pData,A.numRows);
    /****~2.4 X=X+K*(Z-H*X)状态向量后验估计值(结束)****/
  }
  /****~2.使用序贯法进行后验估计结束****/

  PostProcess();//对滤波后的数据进行后处理
}
/**
  ******************************************************************************
  *   函 数 名: Init()
  *   功能说明: 调用A、H、P、Q、R、X的初始化函数对滤波器进行初始化
  *   形    参：
  *   返 回 值: 无
  *   备 注：
  ********************************************************************************   
*/
void CKalmanFilter::Init(const arm_matrix_instance_f32* v_iPara)//对卡尔曼滤波器进行初始化
{
  Init_Matrix();//初始化矩阵
  Init_A();//初始化状态转移A矩阵
  Init_H();//初始化测量转换H矩阵
  Init_P();//初始化状态向量方差P矩阵
  Init_Q();//初始化过程噪声方差Q矩阵
  Init_R();//初始化测量噪声方差R矩阵（假设测量向量相互独立，即R为对角阵）
  Init_X(v_iPara);//初始化状态向量X
}
/*
  以下为虚函数，在子函数中实现
*/
void CKalmanFilter::Init_Matrix(){}//初始化矩阵尺寸
void CKalmanFilter::PreProcess(const arm_matrix_instance_f32* v_iPara,const float32_t* deltaT,const uint8_t* bAvalible){}//执行滤波前预处理
void CKalmanFilter::PostProcess(){}//执行滤波后处理
void CKalmanFilter::Init_A(){}//初始化状态转移A矩阵
void CKalmanFilter::Update_A(){}//更新状态转移A矩阵
void CKalmanFilter::Init_H(){}//初始化测量转换H矩阵
void CKalmanFilter::Update_H(){}//更新测量转换H矩阵
void CKalmanFilter::Init_P(){}//初始化状态向量方差P矩阵
void CKalmanFilter::Init_Q(){}//初始化过程噪声方差Q矩阵
void CKalmanFilter::Update_Q(){}//更新过程噪声方差Q矩阵
void CKalmanFilter::Init_R(){}//初始化测量噪声方差R矩阵（假设测量向量相互独立，即R为对角阵）
void CKalmanFilter::Update_R(){}//更新测量噪声方差R矩阵（假设测量向量相互独立，即R为对角阵）
void CKalmanFilter::Init_X(const arm_matrix_instance_f32* v_iPara){}//初始化状态向量X
void CKalmanFilter::StateEquation(){}//状态转移方程
void CKalmanFilter::MeasuringEquation(){}//测量预测方程
