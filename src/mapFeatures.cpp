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

#include "mapFeatures.h"
#include <Eigen/Dense>
#include <Eigen/Core>
namespace PLSLAM{

// Point features

MapPoint::MapPoint(int idx_, Vector3d point3D_, Mat desc_, int kf_obs_, Vector2d obs_, Vector3d dir_, double sigma2_ ) :
    idx(idx_), point3D(point3D_), inlier(true)
{
    desc_list.push_back( desc_ );
    obs_list.push_back( obs_ );
    kf_obs_list.push_back( kf_obs_ );
    dir_list.push_back( dir_ );
    sigma_list.push_back(sigma2_);
    med_obs_dir = dir_;
    med_desc    = desc_;
}

void MapPoint::addMapPointObservation(Mat desc_, int kf_obs_, Vector2d obs_, Vector3d dir_,  double sigma2_ )
{
    desc_list.push_back( desc_ );
    obs_list.push_back( obs_ );
    kf_obs_list.push_back( kf_obs_ );
    dir_list.push_back( dir_ );
    sigma_list.push_back(sigma2_);
    updateAverageDescDir();
}

void MapPoint::updateAverageDescDir()
{

    // descriptor
    // - check distances between all the observed descriptors
    int n = desc_list.size();
    MatrixXf conf_desc(n,n);
    for(int i = 0; i < n; i++ )
    {
        conf_desc(i,i) = 0;
        for(int j = i+1 ; j < n; j++ )
        {
            int d = norm(desc_list[i],desc_list[j],NORM_HAMMING);
            conf_desc(i,j) = d;
            conf_desc(j,i) = d;
        }
    }

    // - select the one with least mean distance to the rest
    int max_dist = 99999;
    int max_idx  = 0;
    for(int i = 0; i < n; i++)
    {
        vector<int> dist_idx;
        for(int j = 0; j < n; j++)
            dist_idx.push_back( conf_desc(i,j) );
        sort( dist_idx.begin(), dist_idx.end() );
        int idx_median = dist_idx[ int(1+0.5*(n-1)) ];
        if( idx_median < max_dist )
        {
            max_dist = idx_median;
            max_idx  = i;
        }
    }
    med_desc = desc_list[max_idx];

    // direction
    Vector3d med_dir;
    for(int i = 0; i < n; i++)
        med_dir += dir_list[i];
    med_obs_dir = med_dir / n;

}

// Line segment features

MapLine::MapLine(int idx_, Vector6d line3D_, Mat desc_, int kf_obs_, Vector3d obs_, Vector3d dir_, Vector4d pts_, double sigma2_) :
    idx(idx_), line3D(line3D_), inlier(true)
{
    desc_list.push_back( desc_ );
    obs_list.push_back( obs_ );
    kf_obs_list.push_back( kf_obs_ );
    dir_list.push_back( dir_ );
    pts_list.push_back(pts_);
    sigma_list.push_back(sigma2_);
    med_obs_dir = dir_;
    med_desc    = desc_;
}

MapLine::MapLine(int idx_, Vector6d NDw_, Mat desc_, int kf_obs_, Vector4d obs_,  double sigma2) :
    idx(idx_), NDw(NDw_), inlier(true)
{
    desc_list.push_back(desc_);
    NDw_obs_list.push_back(obs_);
    kf_obs_list.push_back(kf_obs_);
    sigma_list.push_back(sigma2);
    med_desc = desc_;
}

void MapLine::addMapLineObservation(Mat desc_, int kf_obs_, Vector3d obs_, Vector3d dir_, Vector4d pts_, double sigma2_)
{
    desc_list.push_back( desc_ );
    obs_list.push_back( obs_ );
    kf_obs_list.push_back( kf_obs_ );
    dir_list.push_back( dir_ );
    pts_list.push_back(pts_);
    sigma_list.push_back(sigma2_);
    updateAverageDescDir();
}

void MapLine::addMapLineObservation(Mat desc_, int kf_obs_, Vector4d obs_, double sigma2_) {
    desc_list.push_back( desc_ );
    NDw_obs_list.push_back( obs_ );
    sigma_list.push_back( sigma2_ );
    kf_obs_list.push_back( kf_obs_ );
    updateAverageDescDir();
}

void MapLine::updateAverageDescDir()
{

    // descriptor
    // - check distances between all the observed descriptors
    int n = desc_list.size();
    MatrixXf conf_desc(n,n);
    for(int i = 0; i < n; i++ )
    {
        conf_desc(i,i) = 0;
        for(int j = i+1 ; j < n; j++ )
        {
            int d = norm(desc_list[i],desc_list[j],NORM_HAMMING);
            conf_desc(i,j) = d;
            conf_desc(j,i) = d;
        }
    }

    // - select the one with least mean distance to the rest
    int max_dist = 99999;
    int max_idx  = 0;
    for(int i = 0; i < n; i++)
    {
        vector<int> dist_idx;
        for(int j = 0; j < n; j++)
            dist_idx.push_back( conf_desc(i,j) );
        sort( dist_idx.begin(), dist_idx.end() );
        int idx_median = dist_idx[ int(1+0.5*(n-1)) ];
        if( idx_median < max_dist )
        {
            max_dist = idx_median;
            max_idx  = i;
        }
    }
    med_desc = desc_list[max_idx];
#ifdef USE_LINE_PLUKER

#else
    // direction
    Vector3d med_dir;
    for(int i = 0; i < n; i++)
        med_dir += dir_list[i];
    med_obs_dir = med_dir / n;
#endif
}

Vector4d MapLine::changePlukerToOrth(const Vector6d& plukerLine)
{
    Vector4d orthLine;
    Matrix3d R = getOrhtRFromPluker(plukerLine);
    Vector3d u1 = R.col(0);
    Vector3d u2 = R.col(1);
    Vector3d u3 = R.col(2);
    orthLine[0] = atan2( u2(2),u3(2) );
    orthLine[1] = asin( -u1(2) );
    orthLine[2] = atan2( u1(1),u1(0) );

    Matrix2d W = getOrthWFromPluker(plukerLine);
    orthLine[3] = asin(W(1,0));

    return orthLine;
}

Vector6d MapLine::changeOrthToPluker(const Vector4d& orthLine)
{
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

Matrix3d MapLine::getOrhtRFromPluker(const Vector6d &plukerLine)
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

Matrix2d MapLine::getOrthWFromPluker(const Vector6d& plukerLine)
{
    Matrix2d temp;
    double nnorm = plukerLine.head(3).norm();
    double dnorm = plukerLine.tail(3).norm();
    double fenmu = sqrt(nnorm*nnorm + dnorm*dnorm);
    temp << nnorm / fenmu, -dnorm / fenmu, dnorm / fenmu, nnorm / fenmu;
    return temp;
}

Matrix<double,6,4> MapLine::jacobianFromPlukerToOrth(const Matrix3d& u, const Matrix2d& w){
    double w1 = w(0,0);
    double w2 = w(1,0);
    Vector3d u1 = u.col(0);
    Vector3d u2 = u.col(1);
    Vector3d u3 = u.col(2);
    Matrix<double,6,4> temp;
    temp.setZero();
    temp.block<3,1>(0,1) = -w1 * u3;
    temp.block<3,1>(0,2) = -w1 * u2;
    temp.block<3,1>(0,3) = -w2 * u1;
    temp.block<3,1>(3,0) = w2 * u3;
    temp.block<3,1>(3,2) = -w2 * u1;
    temp.block<3,1>(3,3) = w1 * u2;
    return temp;
}

}
