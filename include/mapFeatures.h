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
#include <iostream>
#include <vector>
#include <list>
#include <opencv/cv.h>
#include <eigen3/Eigen/Core>

using namespace cv;
using namespace std;
using namespace Eigen;

typedef Matrix<int,5,1>    Vector5i;
typedef Matrix<double,6,1> Vector6d;
typedef Matrix<double,6,6> Matrix6d;

namespace PLSLAM{

class MapPoint
{

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    MapPoint() { }
    MapPoint(int idx_, Vector3d point3D_, Mat desc_, int kf_obs_, Vector2d obs_, Vector3d dir_, double sigma2_ = 1.f);
    ~MapPoint() { }

    void addMapPointObservation(Mat desc_, int kf_obs_, Vector2d obs_, Vector3d dir_,  double sigma2_ = 1.f);
    void updateAverageDescDir();

    int            idx;

    bool           inlier;
    bool           local;
    Vector3d       point3D;
    Vector3d       med_obs_dir;
    Mat            med_desc;

    vector<Mat>      desc_list;       // list with the descriptor of each observation
    vector<Vector2d> obs_list;        // list with the coordinates of each observation
    vector<Vector3d> dir_list;        // list with the direction unit vector of each observation
    vector<int>      kf_obs_list;     // list with KF index where the feature has been observed
    vector<double>   sigma_list;      // list with the sigma scale of each observation

};

class MapLine
{

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    MapLine() { }
    MapLine(int idx_, Vector6d line3D_, Mat desc_, int kf_obs_, Vector3d obs_, Vector3d dir_, Vector4d pts_,  double sigma2_ = 1.f);
    //pluker
    MapLine(int idx_, Vector6d NDw_, Mat desc_, int kf_obs_, Vector4d obs_,  double sigma2_ = 1.f);
    void addMapLineObservation(Mat desc_, int kf_obs_, Vector4d obs_, double sigma2_ = 1.f);


    ~MapLine() { }

    void addMapLineObservation(Mat desc_, int kf_obs_, Vector3d obs_, Vector3d dir_, Vector4d pts_,  double sigma2_ = 1.f);
    void updateAverageDescDir();

    int            idx;

    bool           inlier;
    bool           local;
    Vector6d       line3D;            // 3D endpoints of the line segment
    Vector3d       med_obs_dir;
    Mat            med_desc;

    vector<Mat>      desc_list;       // list with the descriptor of each observation
    vector<Vector3d> obs_list;        // list with the coordinates of each observation ( 2D line equation, normalized by sqrt(lx2+ly2) )
    vector<Vector4d> pts_list;        // list with the coordinates of each endpoint (four coordinates)
    vector<Vector3d> dir_list;        // list with the direction unit vector of each observation (middle point)
    vector<int>      kf_obs_list;     // list with KF index where the feature has been observed
    vector<double>   sigma_list;      // list with the sigma scale of each observation

    //pluker
    Vector6d       NDw;               // pluker coord in the world frame
    Vector4d       orthNDw;           //orth coord for optimize
    vector<Vector4d> NDw_obs_list;
    static Vector4d changePlukerToOrth(const Vector6d& plukerLine);
    static Vector6d changeOrthToPluker(const Vector4d& orthLine );
    static Matrix3d getOrhtRFromPluker(const Vector6d& plukerLine);
    static Matrix2d getOrthWFromPluker(const Vector6d& plukerLine);
    static Matrix<double,6,4> jacobianFromPlukerToOrth(const Matrix3d& u, const Matrix2d& w);

    //debug
    int first_kf_id;
    Matrix4d first_kf_pose;
    Vector4d first_kf_obs;
    double error;
    Vector6d first_NDw;


};

}
