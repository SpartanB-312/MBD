#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>

#define _USE_MATH_DEFINES // 使用math.h中的M_PI宏定义需要

using Eigen::MatrixXd;
using namespace std;

int main()
{
  double g=-9.8;
  double h = 0.001;
  double t = 5;
  int n=(int)t/h;
  double pi=M_PI;
  MatrixXd q(1,n);
  MatrixXd v(1,n);
  MatrixXd a(1,n);

  q(0,0)=pi/3;
  v(0,0)=0;
  a(0,0)=g*sin(q(0,0));
  for(int i=1;i<n;i++)
  {
    q(0,i)=q(0,i-1)+v(0,i-1)*h;
    v(0,i)=v(0,i-1)+a(0,i-1)*h;
    a(0,i)=g*sin(q(0,i));
  }

  ofstream fout("matrixTest.csv");
  for(int i=1;i<n;i++)
  {
    fout << q(0,i) << endl;
  }
	fout.flush();
}