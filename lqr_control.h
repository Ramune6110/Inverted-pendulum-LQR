#ifndef LQR_CONTROL
#define LQR_CONTROL

#include <Eigen/Dense>
using namespace Eigen;

class LQR_control 
{
    private:
        // paramater
        double M;
        double m;
        double J;
        double l;
        double g;
        double alpha;
        // Time 
        double ts;
        double dt;
        double tf;
        // system matrix
        Matrix<double, 4, 4> A;
        Matrix<double, 4, 1> B;
        Matrix<double, 1, 4> C;
        // weight
        Matrix<double, 4, 4> Q;
        Matrix<double, 1, 1> R;
        // Init
        Matrix<double, 4, 1> x;
        Matrix<double, 1, 1> u;
        // Gain
        Matrix<double, 1, 4> K;
    public:
        LQR_control();
        ~LQR_control();
        void simulation();
        MatrixXd calcGainK();
        MatrixXd Runge_Kutta(MatrixXd x);
        MatrixXd model(MatrixXd x, MatrixXd u);
        MatrixXd care(const MatrixXd &A, const MatrixXd &B, const MatrixXd &Q, const MatrixXd &R);
};

#endif // LQR_CONTROL




