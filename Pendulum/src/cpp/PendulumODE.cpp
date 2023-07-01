#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#define _USE_MATH_DEFINES // 使用math.h中的M_PI宏定义需要

using Eigen::MatrixXd;
 
int main()
{
  float g=-9.8;
  float h = 0.001;
  float t = 5;
  int n=(int)t/h;
  float pi=M_PI;
  MatrixXd q(1,n);
  MatrixXd v(1,n);
  MatrixXd a(1,n);

  q(1,1)=pi;
  v(1,1)=0;
  a(1,1)=g*sin(q(1,1));

  for(int i=2;i<=n;i++)
  {
    q(1,i)=q(1,i-1)+v(1,i-1)*h;
    v(1,i)=v(1,i-1)+a(1,i-1)*h;
    a(1,i)=g*sin(q(1,i));
  }
  std::cout << q << std::endl;
}