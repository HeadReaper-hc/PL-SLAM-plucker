/*****************************************************************************
**      Stereo VO and SLAM by combining point and line segment features     **
******************************************************************************
**                                                                          **
**  Copyright(c) 2016-2018, Ruben Gomez-Ojeda, University of Malaga         **
**  Copyright(c) 2016-2018, David Zuñiga-Noël, University of Malaga         **
**  Copyright(c) 2016-2018, MAPIR group, University of Malaga               **
**                                                                          **
**  This program is free software: you can redistribute it and/or modify    **
**  it under the terms of the GNU General Public License (version 3) as     **
**  published by the Free Software Foundation.                              **
**                                                                          **
**  This program is distributed in the hope that it will be useful, but     **
**  WITHOUT ANY WARRANTY; without even the implied warranty of              **
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            **
**  GNU General Public License for more details.                            **
**                                                                          **
**  You should have received a copy of the GNU General Public License       **
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.   **
**                                                                          **
*****************************************************************************/

#pragma once
#include <mutex>
#include <list>
#include <map>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/CholmodSupport>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/Jacobi>
#include <Eigen/src/Core/MatrixBase.h>

#include <g2o/types/slam3d/vertex_se3.h>
#include <g2o/types/slam3d/edge_se3.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/solver.h>
#include <g2o/core/robust_kernel.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/solvers/structure_only/structure_only_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/types/sba/types_six_dof_expmap.h>

#include <DBoW2/TemplatedVocabulary.h>
#include <DBoW2/FORB.h>
#include <DBoW2/BowVector.h>
#include <DBoW2/FClass.h>
#include <DBoW2/FeatureVector.h>
#include <DBoW2/ScoringObject.h>

#include <slamConfig.h>
#include <stereoFrame.h>
#include <stereoFrameHandler.h>
#include <keyFrame.h>
#include <mapFeatures.h>

using namespace std;
using namespace Eigen;
using namespace DBoW2;

typedef Matrix<int,  6,1> Vector6i;
typedef Matrix<float,6,1> Vector6f;
typedef Matrix<float,6,6> Matrix6f;
typedef Matrix<float,7,1> Vector7f;
typedef DBoW2::TemplatedVocabulary<DBoW2::FORB::TDescriptor, DBoW2::FORB> Vocabulary;

namespace PLSLAM
{

class MapHandler
{

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    MapHandler(PinholeStereoCamera* cam_);
    ~MapHandler() { }

    void initialize(KeyFrame* kf0);
    void finishSLAM();
    void addKeyFrame(KeyFrame *curr_kf);

    void addKeyFrame_multiThread(KeyFrame *curr_kf, KeyFrame *prev_kf);
    void handlerThread();

    void startThreads();
    void killThreads();

    void loopClosureThread();
    void localMappingThread();

    int matchKF2KFPoints(KeyFrame *prev_kf, KeyFrame *curr_kf);
    int matchMap2KFPoints();

    int matchKF2KFLines(KeyFrame *prev_kf, KeyFrame *curr_kf);
    int matchMap2KFLines();

    void lookForCommonMatches(KeyFrame *kf0, KeyFrame *&kf1);

    void expandGraphs();
    void formLocalMap();
    void formLocalMap( KeyFrame * kf );
    void formLocalMap_old();
    void removeBadMapLandmarks();
    void removeRedundantKFs();
    void loopClosure();
    bool lookForLoopCandidates(int kf_idx_curr, int &kf_idx_prev);
    void insertKFBowVectorP(KeyFrame *kf);
    void insertKFBowVectorL(KeyFrame *kf);
    void insertKFBowVectorPL(KeyFrame *kf);
    bool isLoopClosure(const KeyFrame* kf0, const KeyFrame* kf1, Vector6d &pose_inc,
                       vector<Vector4i> &lc_pt_idx, vector<Vector4i> &lc_ls_idx,
                       vector<PointFeature*> &lc_points, vector<LineFeature*>  &lc_lines);
    bool computeRelativePoseGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                Vector6d &pose_inc ) const;
    bool computeRelativePoseRobustGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                Vector6d &pose_inc ) const;
    bool loopClosureOptimizationEssGraphG2O();
    bool loopClosureOptimizationCovGraphG2O();
    void loopClosureFuseLandmarks();

    int localBundleAdjustment();
    int levMarquardtOptimizationLBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  );

    //pluker
    int localBundleAdjustmentForPluker();
    int levMarquardtOptimizationLBAForPluker( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  );
    void removeBadMapLandmarksForPluker();
    void localBundleAdjustmentForPlukerWithG2O();

    void globalBundleAdjustment();
    void levMarquardtOptimizationGBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  );

    PinholeStereoCamera* cam;

    vector<KeyFrame*> map_keyframes;
    vector<MapPoint*> map_points;
    vector<MapLine*>  map_lines;

    list<PointFeature*> matched_pt;
    list<LineFeature*>  matched_ls;

    map<int,vector<int>> map_points_kf_idx; // base KF list from which the LM is observed
    map<int,vector<int>> map_lines_kf_idx;

    vector< vector<unsigned int> > full_graph;

    vector< vector<float> > conf_matrix;
    Vocabulary              dbow_voc_p, dbow_voc_l;

    unsigned int max_pt_idx, max_ls_idx, max_kf_idx ;

    KeyFrame *prev_kf, *curr_kf;
    Matrix4d Twf, DT;

    // experiment variables
    Vector7f time;

    // VO status
    mutex m_insert_kf;
    enum VOStatus{
        VO_PROCESSING,
        VO_INSERTING_KF
    };
    VOStatus vo_status;

    // status of the LBA thread
    vector<int> lba_kfs;
    enum LBAState{
        LBA_IDLE,
        LBA_ACTIVE,
        LBA_READY,
        LBA_TERMINATED
    };
    LBAState lba_thread_status;

    // Local Mapping
    std::mutex lba_mutex;
    std::condition_variable lba_start, lba_join;

    vector< Vector3i > lc_idxs,  lc_idx_list;
    vector< Vector6d > lc_poses, lc_pose_list;
    vector< vector<Vector4i> > lc_pt_idxs;
    vector< vector<Vector4i> > lc_ls_idxs;

    std::mutex lc_mutex;
    std::condition_variable lc_start, lc_join;

    enum LCState{
        LC_IDLE,
        LC_ACTIVE,
        LC_READY,
        LC_TERMINATED
    };
    LCState lc_state, lc_thread_status;



    // KF queue
    std::list<pair<KeyFrame*,KeyFrame*>> kf_queue;  // list of curr_kf_mt and prev_kf_mt
    std::mutex kf_queue_mutex;
    std::condition_variable new_kf;
    KeyFrame* curr_kf_mt;
    KeyFrame* prev_kf_mt;

    std::mutex cout_mutex;
    void print_msg(const std::string &msg);

    void SaveKeyFrameTrajectoryTUM(const string &filename);

    static bool lId(KeyFrame* pKF1, KeyFrame* pKF2){
        return pKF1->kf_idx<pKF2->kf_idx;
    }

private:

    bool threads_started;

    inline Matrix3d vectorHat(const Vector3d& vec){
        Matrix3d temp;
        temp << 0, -vec[2], vec[1],
                vec[2], 0, -vec[0],
                -vec[1], vec[0], 0;
        return temp;
    }

    inline Vector6d TransformForPluker(const Matrix4d& T, const Vector6d& vec){
        Matrix6d temp;
        temp.setZero();
        temp.block<3,3>(0,0) = T.block<3,3>(0,0);
        temp.block<3,3>(0,3) = vectorHat(T.block<3,1>(0,3)) * T.block<3,3>(0,0);
        temp.block<3,3>(3,3) = T.block<3,3>(0,0);

        return temp * vec;
    }

    inline Matrix6d getTransformMatrixForPluker(const Matrix4d& T){
        Matrix6d temp;
        temp.setZero();
        temp.block<3,3>(0,0) = T.block<3,3>(0,0);
        temp.block<3,3>(0,3) = vectorHat(T.block<3,1>(0,3)) * T.block<3,3>(0,0);
        temp.block<3,3>(3,3) = T.block<3,3>(0,0);

        return temp;
    }

    inline void updateOrthCoord(Vector4d& D, Vector4d& deltaD, Vector4d& plusD){
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

}
