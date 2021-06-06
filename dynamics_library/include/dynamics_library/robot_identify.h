#ifndef ROBOT_IDENTIFY_H_
#define ROBOT_IDENTIFY_H_

#include <dynamics_library/robot_model.h>
#include <Eigen/Dense>

using namespace Eigen;

namespace robot_dyn
{
class Fourier
{
public:
    Fourier(robot_dyn::RobotModel ROBOT);
    ~Fourier();

    // Which can be redefined by user.
    double wf; 
    double Tf; 
    double h;   //sampling period of system
    // End Which

    unsigned int N;      // order of fourier trajectory
    unsigned int num;    //parameters of each joint

    robot_dyn::RobotModel robot; 

    VectorXd q;
    VectorXd qDot;
    VectorXd qDDot;

private:

};

MatrixXd calcu_Ys(robot_dyn::RobotModel *robot, 
    const VectorXd q, const VectorXd qDot, const VectorXd qDDot);

void calcu_base_dynparam(robot_dyn::RobotModel *robot);

void qr_decompose(robot_dyn::RobotModel *robot, MatrixXd W, MatrixXd &Wb);

// ai_1 ... ai_5, bi_1 ... bi_5
void generate_fourier_trajectory(const VectorXd x, const double t, Fourier *fourier);

double optimal_object_fun(unsigned n, const double *x, double *grad, void *f_data);

void inequality_constraint(unsigned m, double *result, unsigned n, 
    const double *x, double *grad, void *f_data);

void equality_constraint(unsigned m, double *result, unsigned n, 
    const double *x, double *grad, void *f_data);

}

#endif