#ifndef ROBOT_MODEL_H_
#define ROBOT_MODEL_H_

#include <Eigen/Dense>

using namespace Eigen;

namespace robot_dyn
{
Matrix3d Rot(const double alpha, const double theta);

class RobotModel
{
public:
    RobotModel(){};
    ~RobotModel(){};

    void InitModel(const unsigned int DOF);

    void SetKinematicsParameters(const MatrixXd param, const Matrix<Vector3d,1,Dynamic> P_);

    // mi rcxi rcyi rczi Ixxi Ixyi Ixzi Iyyi Iyzi Izzi fvi fsi
    void SetDynamicsParameters(const VectorXd param);

    VectorXd calcu_inv_dyn(const VectorXd q, const VectorXd qDot, const VectorXd qDDot);

    int dof;
    double g;
    VectorXd fv;
    VectorXd fs;

    VectorXd qMin; VectorXd qMax; 
    VectorXd qDotMin; VectorXd qDotMax; 
    VectorXd qDDotMin; VectorXd qDDotMax;

    // mi mrcxi mrcyi mrczi Ixxi Ixyi Ixzi Iyyi Iyzi Izzi fvi fsi
    unsigned int Psi_num;
    unsigned int Ps_num;
    VectorXi Ps_flag;
    VectorXd Ps;
    unsigned int Pb_num;
    VectorXd Pb;
    MatrixXd R1;
    MatrixXd R2;
    double qr_threshold;

private:
    // MDH parameters
    VectorXd alpha;    // x rotation
    VectorXd a;        // x translation
    VectorXd theta;    // z rotation
    VectorXd d;        // z translation
    VectorXd offset;

    Matrix3d I33 = Matrix3d::Identity();

    Vector3d Z;

    VectorXd m;

    Matrix<Matrix3d,1,Dynamic> I;
    Matrix<Matrix3d,1,Dynamic> Ic;
    Matrix<Vector3d,1,Dynamic> Pc;

    Matrix<Vector3d,1,Dynamic> P; 
    Matrix<Matrix3d,1,Dynamic> R;
    Matrix<Matrix3d,1,Dynamic> R_T;

    Matrix<Vector3d,1,Dynamic> w;
    Matrix<Vector3d,1,Dynamic> wDot;
    Matrix<Vector3d,1,Dynamic> vDot;

    Matrix<Vector3d,1,Dynamic> vcDot;
    Matrix<Vector3d,1,Dynamic> F;
    Matrix<Vector3d,1,Dynamic> N;
    
    Matrix<Vector3d,1,Dynamic> f;
    Matrix<Vector3d,1,Dynamic> n;

    VectorXd tau;
};

}

#endif
