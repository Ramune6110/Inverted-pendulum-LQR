#include <iostream>
#include <Eigen/Dense>
#include <cstdio>
#include "lqr_control.h"

using namespace std;
using namespace Eigen;

int main() 
{
    LQR_control lqr;

    lqr.simulation();
    
    return 0;
}

