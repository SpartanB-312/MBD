#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>

#define _USE_MATH_DEFINES // 使用math.h中的M_PI宏定义需要

using Eigen::MatrixXd;
using namespace std;

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

  q(0,0)=pi/3;
  v(0,0)=0;
  a(0,0)=g*sin(q(1,1));

  for(int i=1;i<n;i++)
  {
    q(1,i)=q(1,i-1)+v(1,i-1)*h;
    v(1,i)=v(1,i-1)+a(1,i-1)*h;
    a(1,i)=g*sin(q(1,i));
  }
  //cout << q(1,1) << endl;
  ofstream fout("test.txt");
  fout << q;
  fout.close();
}