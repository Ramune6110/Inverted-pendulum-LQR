#include <iostream>
#include <Eigen/Dense>
#include <cstdio>

using namespace std;
using namespace Eigen;

typedef struct {
    double M;
    double m;
    double J;
    double l;
    double g;
    double alpha;
} param;

typedef struct {
    double ts;
    double dt;
    double tf;
} Time;

typedef struct {
    // system matrix
    Matrix<double, 4, 4> A;
    Matrix<double, 4, 1> B;
    Matrix<double, 1, 4> C;
    // weight
    Matrix<double, 4, 4> Q;
    Matrix<double, 1, 1> R;
    // Init
    Matrix<double, 4, 1> x;
    Matrix<double, 4, 1> xEst;
    Matrix<double, 4, 1> xEst_sec;
    Matrix<double, 4, 4> PEst;
    Matrix<double, 1, 1> u;
    // Gain
    Matrix<double, 1, 4> K;
} matrix;

// prototype
MatrixXd calcGainK(matrix m);
MatrixXd care(const MatrixXd &A, const MatrixXd &B, const MatrixXd &Q, const MatrixXd &R);
MatrixXd Runge_Kutta(matrix m, double dt);
MatrixXd model(MatrixXd x, MatrixXd u, matrix m);

int main() 
{   
    // struct 
    param p;
    Time t;
    matrix m;

    // paramater
    p.M = 0.7;
    p.m = 0.12;
    p.J = 0.009;
    p.l = 0.3;
    p.g = 9.8;
    p.alpha = (p.M + p.m) * (p.J + p.m * pow(p.l, 2.0)) - pow(p.m, 2.0) * pow(p.l, 2.0);

    // Time 
    t.ts = 0.0;
    t.dt = 0.1;
    t.tf = 25.0;

    // Matrix
    m.A << 0, 0, 1, 0,
           0, 0, 0, 1,
           0, -(pow(p.m, 2.0) * pow(p.l, 2.0) * p.g) / p.alpha, 0, 0,
           0, ((p.M + p.m) * p.m * p.g * p.l) / p.alpha, 0, 0;
    
    m.B << 0, 0, (p.J * p.m * pow(p.l, 2.0)) / p.alpha, -(p.m * p.l) / p.alpha;
    m.C << 1, 0, 0, 0;
    m.Q << 10, 0, 0, 0,
           0, 10, 0, 0,
           0, 0, 10, 0,
           0, 0, 0, 10;
    m.R << 1;
    m.x << 1, 0, 0, 0;
    m.xEst << 1, 0, 0, 0;
    m.xEst_sec << 1, 0, 0, 0;
    m.PEst << 0.1, 0, 0, 0,
              0, 0.1, 0, 0,
              0, 0, 0.1, 0,
              0, 0, 0, 0.1;
    m.u << 0;

    // save data
    FILE *fp;
    if ((fp = fopen("data.txt", "w")) == NULL) {
        printf("Error\n");
        exit(1);
    }
    fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", t.ts, m.x(0), m.x(1), m.x(2), m.x(3));
    // main loop
    for (t.ts; t.ts <= t.tf; t.ts += t.dt) {
        m.K = calcGainK(m);
        m.u = -m.K * m.x;
        m.x = Runge_Kutta(m, t.dt);
        fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", t.ts, m.x(0), m.x(1), m.x(2), m.x(3));
    }
    fclose(fp);

    return 0;
}

MatrixXd calcGainK(matrix m) {
    // LQR
    MatrixXd  S = care(m.A, m.B, m.Q, m.R);
    return m.R.inverse() * (m.B.transpose() * S);
}

MatrixXd care(const MatrixXd &A, const MatrixXd &B, const MatrixXd &Q, const MatrixXd &R) {
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

MatrixXd Runge_Kutta(matrix m, double dt) {
    // -- runge-kutta --
    MatrixXd k1 = model(m.x, m.u, m);
    MatrixXd k2 = model(m.x + k1 * dt / 2.0, m.u, m);
    MatrixXd k3 = model(m.x + k2 * dt / 2.0, m.u, m);
    MatrixXd k4 = model(m.x + k3 * dt, m.u, m);

    return m.x = m.x + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt / 6.0;
}

MatrixXd model(MatrixXd x, MatrixXd u, matrix m) {
    return x = m.A * x + m.B * u;
}


