#include <iostream>
#include <Eigen/Dense>
#include <cstdio>
#include "lqr_control.h"

using namespace std;
using namespace Eigen;

LQR_control::LQR_control() {
    // paramater
    M = 0.7;
    m = 0.12;
    J = 0.009;
    l = 0.3;
    g = 9.8;
    alpha = (M + m) * (J + m * pow(l, 2.0)) - pow(m, 2.0) * pow(l, 2.0);
    // Time 
    ts = 0.0;
    dt = 0.1; 
    tf = 25.0;
    // Matrix
    A << 0, 0, 1, 0,
         0, 0, 0, 1,
         0, -(pow(m, 2.0) * pow(l, 2.0) * g) / alpha, 0, 0,
         0, ((M + m) * m * g * l) / alpha, 0, 0;

    B << 0, 0, (J + m * pow(l, 2.0)) / alpha, -(m * l) / alpha;

    C << 1, 0, 0, 0;

    Q << 10, 0, 0, 0,
         0, 10, 0, 0,
         0, 0, 10, 0,
         0, 0, 0, 10;

    R << 1;

    x << 1, 0, 0, 0;
}

LQR_control::~LQR_control() {
    cout << "Finish" << endl;
}

MatrixXd LQR_control::calcGainK() {
    // LQR
    MatrixXd  S = care(A, B, Q, R);
    return R.inverse() * (B.transpose() * S);
}

MatrixXd LQR_control::care(const MatrixXd &A, const MatrixXd &B, const MatrixXd &Q, const MatrixXd &R) {
      const size_t dim_x = A.rows();
      
      // Set Hamilton Matrix
      Eigen::MatrixXd Ham(2*dim_x, 2*dim_x);
      Ham << A, -B*R.inverse()*B.transpose(), -Q, -A.transpose();

      // calc eigenvalues and eigenvectors
      Eigen::EigenSolver<Eigen::MatrixXd> Eigs(Ham);
      if (Eigs.info() != Eigen::Success) abort();

      // extract stable eigenvectors into 'eigvec'
      Eigen::MatrixXcd eigvec(2*dim_x, dim_x);
      int j = 0;

      // store those with negative real number
      for(int i = 0; i < 2*dim_x; ++i){
          if(Eigs.eigenvalues()[i].real() < 0){
              eigvec.col(j) = Eigs.eigenvectors().block(0, i, 2*dim_x, 1);
              ++j;
          }
      }

      // calc S with stable eigen vector matrix
      Eigen::MatrixXcd U(dim_x, dim_x);
      Eigen::MatrixXcd V(dim_x, dim_x);

      U = eigvec.block(0,0,dim_x,dim_x);
      V = eigvec.block(dim_x,0,dim_x,dim_x);

      return (V * U.inverse()).real();
}

MatrixXd LQR_control::Runge_Kutta(MatrixXd x) {
    // -- runge-kutta --
    MatrixXd k1 = model(x, u);
    MatrixXd k2 = model(x + k1 * dt / 2.0, u);
    MatrixXd k3 = model(x + k2 * dt / 2.0, u);
    MatrixXd k4 = model(x + k3 * dt, u);

    return x = x + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt / 6.0;
}

MatrixXd LQR_control::model(MatrixXd x, MatrixXd u) {
    return x = A * x + B * u;
}

void LQR_control::simulation() {
    // save data
    FILE *fp;
    if ((fp = fopen("data_class.txt", "w")) == NULL) {
        printf("Error\n");
        exit(1);
    }
    fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", ts, x(0), x(1), x(2), x(3));
    // main loop
    for (ts; ts <= tf; ts += dt) {
        K = calcGainK();
        u = -K * x;
        x = Runge_Kutta(x);
        fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", ts, x(0), x(1), x(2), x(3));
    }
    fclose(fp);
}
