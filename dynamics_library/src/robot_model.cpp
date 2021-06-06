#include <dynamics_library/robot_model.h>
#include <dynamics_library/robot_math.h>
#include <math.h>

namespace robot_dyn
{
Matrix3d Rot(const double alpha, const double theta)
{
    Matrix3d R;

    R << cos(theta),              -sin(theta),                0.0,
         cos(alpha)*sin(theta),    cos(alpha)*cos(theta),    -sin(alpha),
         sin(alpha)*sin(theta),    sin(alpha)*cos(theta),     cos(alpha);

    return R;
}

void RobotModel::InitModel(const unsigned int DOF)
{
    dof = DOF;
    g = -9.81;
    fv.resize(dof);
    fs.resize(dof);

    qMin.resize(dof);
    qMax.resize(dof);
    qDotMin.resize(dof);
    qDotMax.resize(dof);
    qDDotMin.resize(dof);
    qDDotMax.resize(dof);

    Psi_num = 12;
    Ps_num = Psi_num*dof;
    Ps_flag.resize(Ps_num);
    Ps.resize(Ps_num);
    Pb_num = 0;
    qr_threshold = 1.0e-10;

    alpha.resize(dof);
    a.resize(dof);
    theta.resize(dof);
    d.resize(dof);
    offset.resize(dof);

    Z << 0.0, 0.0, 1.0;

    m.resize(dof);

    I.resize(dof);
    Ic.resize(dof);
    Pc.resize(dof);

    P.resize(dof+1);
    R.resize(dof+1);
    R_T.resize(dof+1);

    w.resize(dof+1);
    wDot.resize(dof+1);
    vDot.resize(dof+1);

    vcDot.resize(dof);
    F.resize(dof);
    N.resize(dof);

    f.resize(dof+1);
    n.resize(dof+1);

    tau.resize(dof);
}

void RobotModel::SetKinematicsParameters(const MatrixXd param, const Matrix<Vector3d,1,Dynamic> P_)
{
    alpha = param.row(0).transpose();
    a = param.row(1).transpose();
    d = param.row(2).transpose();
    offset = param.row(3).transpose();

    for (unsigned int j=0; j<dof; j++)
    {
        P(j) = P_(j);
    }

    P(dof) << 0.0, 0.0, 0.0;
}

// mi rcxi rcyi rczi Ixxi Ixyi Ixzi Iyyi Iyzi Izzi fvi fsi
void RobotModel::SetDynamicsParameters(const VectorXd param)
{
    for (unsigned int i=0; i<dof; i++)
    {
        m(i) = param(i*Psi_num+0);

        Pc(i) << param(i*Psi_num+1), param(i*Psi_num+2), param(i*Psi_num+3);
    
        I(i) << param(i*Psi_num+4), param(i*Psi_num+5), param(i*Psi_num+6),
                param(i*Psi_num+5),  param(i*Psi_num+7), param(i*Psi_num+8),
                param(i*Psi_num+6), param(i*Psi_num+8),  param(i*Psi_num+9);

        fv(i) = param(i*Psi_num+10);
        fs(i) = param(i*Psi_num+11);

        Ic(i) = I(i)-m(i)*(Pc(i).transpose()*Pc(i)*I33-Pc(i)*Pc(i).transpose());
    }
}

VectorXd RobotModel::calcu_inv_dyn(const VectorXd q, const VectorXd qDot, const VectorXd qDDot)
{
    theta = q+offset;

    for (unsigned int i=0; i<dof; i++)
    {
        R(i) = Rot(alpha(i), theta(i));
        R_T(i) = R(i).transpose();
    }

    R(dof) = Rot(0.0, 0.0);
    R_T(dof) = R(dof).transpose();

    w(0) << 0.0, 0.0, 0.0;
    wDot(0) << 0.0, 0.0, 0.0;
    vDot(0) << 0.0, 0.0, g;
    
    f(dof) << 0.0, 0.0, 0.0;
    n(dof) << 0.0, 0.0, 0.0;

    for (unsigned int i=1; i<=dof; i++)
    {
        w(i) = R(i-1)*w(i-1)+qDot(i-1)*Z;
        wDot(i) = R(i-1)*wDot(i-1)+R(i-1)*w(i-1).cross(qDot(i-1)*Z)+qDDot(i-1)*Z;
        vDot(i) = R(i-1)*(wDot(i-1).cross(P(i-1))+w(i-1).cross(w(i-1).cross(P(i-1)))+vDot(i-1));
        vcDot(i-1) = wDot(i).cross(Pc(i-1))+w(i).cross(w(i).cross(Pc(i-1)))+vDot(i);
        F(i-1) = m(i-1)*vcDot(i-1);
        N(i-1) = Ic(i-1)*wDot(i)+w(i).cross(Ic(i-1)*w(i));
    }
    
    for (int i=dof-1; i>=0; i--)
    {
        f(i) = R_T(i+1)*f(i+1)+F(i);
        n(i) = N(i)+R_T(i+1)*n(i+1)+Pc(i).cross(F(i))+P(i+1).cross(R_T(i+1)*f(i+1));
        tau(i) = n(i).transpose()*Z+fv(i)*qDot(i)+fs(i)*sign(qDot(i));
    }

    return tau;
}

}









