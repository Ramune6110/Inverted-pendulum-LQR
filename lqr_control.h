#ifndef LQR_CONTROL
#define LQR_CONTROL

#include <Eigen/Dense>
using namespace Eigen;

class LQR_control 
{
    // LQR paramater
    private:
        // system paramater
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
        // state and input
        Matrix<double, 4, 1> x;
        Matrix<double, 1, 1> u;
        // LQR Gain
        Matrix<double, 1, 4> K;
    // Kalmanfilter parametar
    private:
        // noise
        Matrix<double, 1, 1> Qsigma;
        Matrix<double, 1, 1> Rsigma;
        // Observation
        Matrix<double, 1, 1> z;
        // Kalman Gain
        Matrix<double, 4, 1> G;
        // state estimation
        Matrix<double, 1, 1> uEst;
        Matrix<double, 1, 1> zEst;
        Matrix<double, 4, 1> xEst;
        Matrix<double, 4, 4> PEst;
        Matrix<double, 4, 1> xPred;
        Matrix<double, 4, 4> PPred;
        // estimate integral
        Matrix<double, 4, 1> XEST;
    // LQR function
    private:
        MatrixXd calcGainK();
        MatrixXd Runge_Kutta(MatrixXd x, MatrixXd u);
        MatrixXd model(MatrixXd x, MatrixXd u);
        MatrixXd care(const MatrixXd &A, const MatrixXd &B, const MatrixXd &Q, const MatrixXd &R);
    // Kalmanfilter function
    private:
        void Kalmanfilter();
    public:
        LQR_control();
        ~LQR_control();
        void simulation();
};

#endif // LQR_CONTROL




