//
// Created by hc on 2020/4/2.
//

#ifndef PL_SLAM_G2O_TYPES_H
#define PL_SLAM_G2O_TYPES_H

#include <g2o/core/base_vertex.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/base_multi_edge.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/core/block_solver.h>
#include <iostream>
using namespace g2o;
typedef g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType> SlamLinearSolver;

inline Matrix3d vechat(const Vector3d& vec){
    Matrix3d temp;
    temp << 0, -vec[2], vec[1],
            vec[2], 0, -vec[0],
            -vec[1], vec[0], 0;
    return temp;
}


//point vertex
class VertexLMPointXYZ : public BaseVertex<3, Vector3d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VertexLMPointXYZ() : BaseVertex<3,Vector3d>() {};
    virtual bool read(std::istream& is) { return true; };
    virtual bool write(std::ostream& os) const { return true; };

    virtual void setToOriginImpl() {
        _estimate.fill(0.);
    }

    virtual void oplusImpl(const double* update)
    {
        Eigen::Map<const Vector3d> v(update);
        _estimate += v;
    }

    virtual int estimateDimension() const{
        return 3;
    }
};

//line vertex
class VertexLMLineOrth : public BaseVertex<4, Vector4d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VertexLMLineOrth() : BaseVertex<4,Vector4d>() {};
    virtual bool read(std::istream& is) { return true; };
    virtual bool write(std::ostream& os) const { return true; };

    virtual void setToOriginImpl() {
        _estimate.fill(0.);
    }

    virtual void oplusImpl(const double* update)
    {
        Eigen::Map<const Vector4d> delta(update);
        Vector4d plusD;
        updateOrthCoord(_estimate, delta, plusD);
        _estimate = plusD;
    }

    inline void updateOrthCoord(const Vector4d& D, const Vector4d& deltaD, Vector4d& plusD){
        // ref: 2001, Adrien Bartol,Peter Sturm ,Structure-From-Motion Using Lines: Representation, Triangulation and Bundle Adjustment

        // theta --> U,  phi --> W
        Eigen::Vector3d theta = D.head(3);

        double phi = D(3);
        //Vector3d theta = orth.head(3);
        //double phi = orth[3];
        double s1 = sin(theta[0]);
        double c1 = cos(theta[0]);
        double s2 = sin(theta[1]);
        double c2 = cos(theta[1]);
        double s3 = sin(theta[2]);
        double c3 = cos(theta[2]);
        Eigen::Matrix3d R;
        R <<
          c2 * c3,   s1 * s2 * c3 - c1 * s3,   c1 * s2 * c3 + s1 * s3,
                c2 * s3,   s1 * s2 * s3 + c1 * c3,   c1 * s2 * s3 - s1 * c3,
                -s2,                  s1 * c2,                  c1 * c2;
        double w1 = cos(phi);
        double w2 = sin(phi);

        // update
        Eigen::Vector3d _delta_theta = deltaD.head(3);
        double _delta_phi = deltaD(3);
        Eigen::Matrix3d Rz;
        Rz << cos(_delta_theta(2)), -sin(_delta_theta(2)), 0,
                sin(_delta_theta(2)), cos(_delta_theta(2)), 0,
                0, 0, 1;

        Eigen::Matrix3d Ry;
        Ry << cos(_delta_theta(1)), 0., sin(_delta_theta(1)),
                0., 1., 0.,
                -sin(_delta_theta(1)), 0., cos(_delta_theta(1));

        Eigen::Matrix3d Rx;
        Rx << 1., 0., 0.,
                0., cos(_delta_theta(0)), -sin(_delta_theta(0)),
                0., sin(_delta_theta(0)), cos(_delta_theta(0));
        R = R * Rx * Ry * Rz;

        Eigen::Matrix2d W;
        W << w1, -w2, w2, w1;
        Eigen::Matrix2d delta_W;
        delta_W << cos(_delta_phi), -sin(_delta_phi),sin(_delta_phi), cos(_delta_phi);
        W = W * delta_W;

        // U' -- > theta'. W' --> phi'

        Eigen::Vector3d u1 = R.col(0);
        Eigen::Vector3d u2 = R.col(1);
        Eigen::Vector3d u3 = R.col(2);
        plusD[0] = atan2( u2(2),u3(2) );
        plusD[1] = asin( -u1(2) );
        plusD[2] = atan2( u1(1),u1(0) );

        plusD[3] = asin( W(1,0) );

//////////////////////////////////////////////////////////
        /*
        // SO3参数方法，得到的雅克比更上面一样的。用上面的形式就OK。
        Eigen::Map<const Eigen::Vector3d> theta(x);
        double phi = *(x + 3);
        double s1 = sin(theta[0]);
        double c1 = cos(theta[0]);
        double s2 = sin(theta[1]);
        double c2 = cos(theta[1]);
        double s3 = sin(theta[2]);
        double c3 = cos(theta[2]);
        Matrix3d R;
        R <<
          c2 * c3,   s1 * s2 * c3 - c1 * s3,   c1 * s2 * c3 + s1 * s3,
                c2 * s3,   s1 * s2 * s3 + c1 * c3,   c1 * s2 * s3 - s1 * c3,
                -s2,                  s1 * c2,                  c1 * c2;

        Sophus::SO3<double> U = Sophus::SO3<double>::exp(theta);
        Sophus::SO3<double> U1(R);

        std::cout << U.matrix() << "\n\n" <<U1.matrix()<<"\n\n"<<R<<"\n\n";

        std::cout << theta <<"\n\n" << U1.log() << "\n\n"<<  Sophus::SO3<double>::exp(U1.log()).matrix() << "\n\n";
         */
    }
};

//Pose Vertex
class VertexLMPose : public BaseVertex<6, Matrix4d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VertexLMPose() : BaseVertex<6,Matrix4d>() {}

    virtual bool read(std::istream& is) { return true; };
    virtual bool write(std::ostream& os) const { return true; };

    virtual void setToOriginImpl() {
        _estimate.setIdentity();
    }

    virtual void oplusImpl(const double* update)  //左乘 平移在前，旋转在后
    {
        Eigen::Map<const Vector6d> delta(update);

        Vector3d omega = delta.tail(3);

        Vector3d rot = delta.tail(3);
        double theta = rot.norm();
        double half_theta = 0.5*(theta);

        double imag_factor;
        double real_factor = cos(half_theta);
        if((theta)<1e-10)
        {
            double theta_sq = (theta)*(theta);
            double theta_po4 = theta_sq*theta_sq;
            imag_factor = 0.5-0.0208333*theta_sq+0.000260417*theta_po4;
        }
        else
        {
            double sin_half_theta = sin(half_theta);
            imag_factor = sin_half_theta/(theta);
        }

        Quaterniond temp(real_factor,
                         imag_factor*omega.x(),
                         imag_factor*omega.y(),
                         imag_factor*omega.z());
        Matrix3d delta_rot = temp.toRotationMatrix();
        _estimate.block<3,3>(0,0) = delta_rot * _estimate.block<3,3>(0,0);
        _estimate.block<3,1>(0,3) += delta.head(3);
    }
};

class EdgePosePoint : public BaseBinaryEdge<2, Vector2d, VertexLMPointXYZ, VertexLMPose>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgePosePoint() : BaseBinaryEdge<2, Vector2d, VertexLMPointXYZ, VertexLMPose>() {}

    bool read(std::istream &is) { return true; }

    bool write(std::ostream &os) const { return true; }

    void SetParams(const double &fx_, const double &fy_, const double &cx_, const double &cy_) {
        fx = fx_;
        fy = fy_;
        cx = cx_;
        cy = cy_;
    }

    void computeError() {
        Vector3d Pc = computePc();
        Vector2d obs(_measurement);

        _error = obs - cam_project(Pc);
        //std::cout<<"_error: "<<_error<<std::endl;
    }

    Vector3d computePc() {
        const VertexLMPointXYZ *vPoint = dynamic_cast<const VertexLMPointXYZ *>(_vertices[0]);
        const VertexLMPose *vPose = dynamic_cast<const VertexLMPose *>(_vertices[1]);

        const Matrix4d &ns = vPose->estimate();
        Matrix3d Rcw = ns.block<3,3>(0,0);
        Vector3d Pcw = ns.block<3,1>(0,3);
        const Vector3d &Pw = vPoint->estimate();

        Vector3d Pc = Rcw * Pw + Pcw;

        return Pc;
        //Vector3d Pwc = Rwb*Pbc + Pwb;
        //Matrix3d Rcw = (Rwb*Rbc).transpose();
        //Vector3d Pcw = -Rcw*Pwc;
        //Vector3d Pc = Rcw*Pw + Pcw;
    }

    inline Vector2d project2d(const Vector3d &v) const {
        Vector2d res;
        res(0) = v(0) / v(2);
        res(1) = v(1) / v(2);
        return res;
    }

    Vector2d cam_project(const Vector3d &trans_xyz) const {
        Vector2d proj = project2d(trans_xyz);
        Vector2d res;
        res[0] = proj[0] * fx + cx;
        res[1] = proj[1] * fy + cy;
        return res;
    }

    bool isDepthPositive() {
        Vector3d Pc = computePc();
        return Pc(2) > 0.0;
    }

    //
    virtual void linearizeOplus(){
        const VertexLMPointXYZ *vPoint = dynamic_cast<const VertexLMPointXYZ *>(_vertices[0]);
        const VertexLMPose *vPose = dynamic_cast<const VertexLMPose *>(_vertices[1]);

        const Matrix3d Rot = vPose->estimate().block<3,3>(0,0);
        const Vector3d Trans = vPose->estimate().block<3,1>(0,3);

        const Vector3d Pw = vPoint->estimate();

        Vector3d Pc = computePc();

        double x = Pc(0);
        double y = Pc(1);
        double z = Pc(2);
        double invz = 1.0 / z;
        double invz2 = invz * invz;

        Matrix<double,2,3> jac_pixel_cam;
        jac_pixel_cam << fx / z, 0, -fx * x * invz2,
                         0, fy / z, -fy * y * invz2;

        _jacobianOplusXi = - jac_pixel_cam * Rot;

        _jacobianOplusXj.block<2,3>(0,0) = - jac_pixel_cam;
        _jacobianOplusXj.block<2,3>(0,3) = - jac_pixel_cam * (- ::vechat(Rot*Pw));
    }

protected:
    double fx, fy, cx, cy;
};

class EdgePoseLine : public BaseBinaryEdge<4, Vector4d, VertexLMLineOrth, VertexLMPose>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgePoseLine() : BaseBinaryEdge<4, Vector4d, VertexLMLineOrth, VertexLMPose>() {}

    bool read(std::istream &is) { return true; }

    bool write(std::ostream &os) const { return true; }

    void SetParams(const double &fx_, const double &fy_, const double &cx_, const double &cy_) {
        fx = fx_;
        fy = fy_;
        cx = cx_;
        cy = cy_;
    }

    void computeError() {
        const VertexLMLineOrth* vLine = dynamic_cast<const VertexLMLineOrth* >(_vertices[0]);
        const VertexLMPose* vPose = dynamic_cast<const VertexLMPose* >(_vertices[1]);

        const Vector4d obs(_measurement);

        const Matrix4d Tcw = vPose->estimate();
        const Vector4d Lw = vLine->estimate();

        Vector6d plukerLw = changeOrthToPluker(Lw);
        Vector6d plukerLc = getTransformMatrix(Tcw) * plukerLw;

        Vector3d plukerLcPixel = getPlukerK() * plukerLc.head(3);
        double lx = plukerLcPixel(0);
        double ly = plukerLcPixel(1);
        double lz = plukerLcPixel(2);

        double fenmu = sqrt( lx*lx + ly*ly);

        Vector4d error;
        error(0) = ( lx * obs(0) + ly * obs(1) + lz ) / fenmu;
        error(1) = ( lx * obs(2) + ly * obs(3) + lz ) / fenmu;
        error(2) = 0;
        error(3) = 0;

        _error = error;
       // std::cout<<"_error: "<<_error.transpose()<<std::endl;
    }

    Matrix3d getPlukerK(){
        Matrix3d plukerK;
        plukerK << fy,     0,      0,
                0,      fx,     0,
                -fy*cx, -fx*cy, fx*fy;
        return plukerK;
    }

    Matrix6d getTransformMatrix(const Matrix4d& T){
        Matrix6d temp;
        temp.setZero();
        temp.block<3,3>(0,0) = T.block<3,3>(0,0);
        temp.block<3,3>(0,3) = ::vechat(T.block<3,1>(0,3)) * T.block<3,3>(0,0);
        temp.block<3,3>(3,3) = T.block<3,3>(0,0);

        return temp;
    }

    Vector6d changeOrthToPluker(const Vector4d& orthLine ){
        double s1 = sin(orthLine[0]);
        double c1 = cos(orthLine[0]);
        double s2 = sin(orthLine[1]);
        double c2 = cos(orthLine[1]);
        double s3 = sin(orthLine[2]);
        double c3 = cos(orthLine[2]);
        Eigen::Matrix3d R;
        R <<
          c2 * c3,   s1 * s2 * c3 - c1 * s3,   c1 * s2 * c3 + s1 * s3,
                c2 * s3,   s1 * s2 * s3 + c1 * c3,   c1 * s2 * s3 - s1 * c3,
                -s2,                  s1 * c2,                  c1 * c2;

        double w1 = cos(orthLine[3]);
        double w2 = sin(orthLine[3]);

        Vector6d plukerLine;
        plukerLine.head(3) = w1 * R.col(0);
        plukerLine.tail(3) = w2 * R.col(1);
        return plukerLine;
    }

    virtual void linearizeOplus() {
        const VertexLMLineOrth* vLine = dynamic_cast<const VertexLMLineOrth* >(_vertices[0]);
        const VertexLMPose* vPose = dynamic_cast<const VertexLMPose* >(_vertices[1]);

        const Vector4d obs(_measurement);

        const Matrix4d Tcw = vPose->estimate();
        const Vector4d Lw = vLine->estimate();

        const Matrix3d Rcw = Tcw.block<3,3>(0,0);
        const Vector3d Pcw = Tcw.block<3,1>(0,3);

        Vector6d plukerLw = changeOrthToPluker(Lw);
        Vector6d plukerLc = getTransformMatrix(Tcw) * plukerLw;

        Vector3d plukerLcPixel = getPlukerK() * plukerLc.head(3);
        double lx = plukerLcPixel(0);
        double ly = plukerLcPixel(1);
        double lz = plukerLcPixel(2);
        double fenmu = sqrt(lx*lx + ly*ly);

        Vector2d error;
        error(0) = ( lx * obs(0) + ly * obs(1) + lz ) / fenmu;
        error(1) = ( lx * obs(2) + ly * obs(3) + lz ) / fenmu;

        Matrix<double,1,3> jac_err0_lcPixel;
        jac_err0_lcPixel << -lx * error(0) / (fenmu*fenmu) + obs(0) /fenmu,
                            -ly * error(0) / (fenmu*fenmu) + obs(1) /fenmu,
                            1.0 / fenmu;
        Matrix<double,1,3> jac_err1_lcPixel;
        jac_err1_lcPixel << -lx * error(1) / (fenmu*fenmu) + obs(2) /fenmu,
                -ly * error(1) / (fenmu*fenmu) + obs(3) /fenmu,
                1.0 / fenmu;

        Matrix<double,3,6> jac_lcPixel_lc;
        jac_lcPixel_lc.setZero();
        jac_lcPixel_lc.block<3,3>(0,0) = getPlukerK();

        Matrix<double,6,6> jac_lc_rt;
        jac_lc_rt.setZero();
        jac_lc_rt.block<3,3>(0,0) = - vechat(Rcw*Lw.tail(3));
        jac_lc_rt.block<3,3>(0,3) = - vechat(Rcw*Lw.head(3)) - vechat(Pcw) * vechat(Rcw*Lw.tail(3));

        Matrix<double,4,6> jac_error_pose;
        jac_error_pose.block<1,6>(0,0) = jac_err0_lcPixel * jac_lcPixel_lc * jac_lc_rt;
        jac_error_pose.block<1,6>(1,0) = jac_err1_lcPixel * jac_lcPixel_lc * jac_lc_rt;
        jac_error_pose.block<2,6>(2,0) = Matrix<double,2,6>::Zero();
        _jacobianOplusXj = jac_error_pose;

        Matrix<double,6,6> jac_lc_lw;
        jac_lc_lw = getTransformMatrix(Tcw);

        Matrix3d U = getOrhtRFromPluker(plukerLw);
        Matrix2d W = getOrthWFromPluker(plukerLw);

        Matrix<double,6,4> jac_lw_orth;
        jac_lw_orth = jacobianFromPlukerToOrth(U, W);

        Matrix<double,4,4> jac_error_orth;
        jac_error_orth.block<1,4>(0,0) = jac_err0_lcPixel * jac_lcPixel_lc * jac_lc_lw * jac_lw_orth;
        jac_error_orth.block<1,4>(1,0) = jac_err1_lcPixel * jac_lcPixel_lc * jac_lc_lw * jac_lw_orth;
        jac_error_orth.block<2,4>(2,0) = Matrix<double,2,4>::Zero();
        _jacobianOplusXi = jac_error_orth;

    }

    Matrix<double,6,4> jacobianFromPlukerToOrth(const Matrix3d& u, const Matrix2d& w){
        double w1 = w(0,0);
        double w2 = w(1,0);
        Vector3d u1 = u.col(0);
        Vector3d u2 = u.col(1);
        Vector3d u3 = u.col(2);
        Matrix<double,6,4> temp;
        temp.setZero();
        temp.block<3,1>(0,1) = -w1 * u3;
        temp.block<3,1>(0,2) = w1 * u2;
        temp.block<3,1>(0,3) = -w2 * u1;
        temp.block<3,1>(3,0) = w2 * u3;
        temp.block<3,1>(3,2) = -w2 * u1;
        temp.block<3,1>(3,3) = w1 * u2;
        return temp;
    }

    Matrix2d getOrthWFromPluker(const Vector6d& plukerLine)
    {
        Matrix2d temp;
        double nnorm = plukerLine.head(3).norm();
        double dnorm = plukerLine.tail(3).norm();
        double fenmu = sqrt(nnorm*nnorm + dnorm*dnorm);
        temp << nnorm / fenmu, -dnorm / fenmu, dnorm / fenmu, nnorm / fenmu;
        return temp;
    }

    Matrix3d getOrhtRFromPluker(const Vector6d &plukerLine)
    {
        Matrix3d R;
        Vector3d n = plukerLine.head(3);
        Vector3d d = plukerLine.tail(3);
        Vector3d n0 = n;
        Vector3d d0 = d;
        n.normalize();
        d.normalize();
        R.col(0) = n;
        R.col(1) = d;
        R.col(2) = n0.cross(d0) / (n0.cross(d0).norm());
        return R;
    }



protected:
    double fx,fy,cx,cy;

};



#endif //PL_SLAM_G2O_TYPES_H
