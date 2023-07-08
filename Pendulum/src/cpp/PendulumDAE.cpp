#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>
#include <ctime>

#define _USE_MATH_DEFINES // 使用math.h中的M_PI宏定义需要

using Eigen::MatrixXd;
using namespace std;

int main()
{
  clock_t start = clock();

  double g=-9.8;
  double h = 0.001;
  double t = 5;
  int n=(int)t/h;
  double pi=M_PI;
  double alpha=25;
  double beta=50;

  MatrixXd q(3,n);
  MatrixXd v(3,n);
  MatrixXd a(3,n);
  MatrixXd M(3,3);//mass
  M << 1, 0, 0,
       0, 1, 0,
       0, 0, 0;
  MatrixXd B(3,1);//force
  B << 0,
       M(1,1)*g,
       0;
  MatrixXd Phi(2,1);
  MatrixXd Phiq(2,3);
  MatrixXd PhiT(2,1);
  MatrixXd P(2,1);
  MatrixXd P1(2,1);//for Baumgrate
  MatrixXd LHS(5,5);
  MatrixXd RHS(5,1);
  MatrixXd zeros4LHS = MatrixXd::Zero(2,2);
  MatrixXd sol(5,1);

  q(0,0)=0.5;q(1,0)=-0.5*sqrt(3);q(2,0)=pi/6;
  v.col(0) << 0,
              0,
              0;
  for(int i=0;i<n-1;i++)
  {
    P << - v(2,i)*v(2,i)*cos(q(2,i)),
           v(2,i)*v(2,i)*sin(q(2,i));
    Phi << q(0,i)-sin(q(2,i)),
           q(1,i)+cos(q(2,i));
    Phiq << 1, 0, -cos(q(2,i)),
            0, 1, -sin(q(2,i));
    PhiT=Phiq*(v.col(i));
    P1=P-2*alpha*PhiT-(beta*beta)*Phi;
    LHS << M, Phiq.transpose(),
           Phiq, zeros4LHS;
    RHS << B,
           P1;
    sol = LHS.inverse()*RHS;
    a(0,i)=sol(0,0);a(1,i)=sol(1,0);a(2,i)=sol(2,0);
    q.col(i+1)=q.col(i)+v.col(i)*h;
    v.col(i+1)=v.col(i)+a.col(i)*h;
    //q(0,i+1)=q(0,i)+v(0,i)*h;
    //v(0,i+1)=v(0,i)+a(0,i)*h;
  }
  //output
  ofstream fout("PendulumDAE.csv");
  for(int i=0;i<n-1;i++)
  {
    fout << q(2,i) << endl;
  }
	fout.flush();

  clock_t end = clock();
  cout << (double)(end - start) / CLOCKS_PER_SEC << 's' << endl;
  getchar();
  getchar();
  return 0;
}