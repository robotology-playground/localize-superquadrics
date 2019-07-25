/******************************************************************************
 *                                                                            *
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file nlp.cpp
 * @authors: Ugo Pattacini <ugo.pattacini@iit.it>
 */

#include <cmath>
#include <limits>

#include <yarp/math/Math.h>
#include <yarp/math/SVD.h>

#include "nlp_vertical.h"

using namespace std;
using namespace yarp::sig;
using namespace yarp::math;

/****************************************************************/
bool SuperQuadricVerticalNLP::get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
                                   Ipopt::Index &nnz_jac_g,
                                   Ipopt::Index &nnz_h_lag,
                                   IndexStyleEnum &index_style)
{
    n=6;
    m = nnz_jac_g = nnz_h_lag = 0;
    index_style=TNLP::C_STYLE;
    return true;
}

/****************************************************************/
bool SuperQuadricVerticalNLP::get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l,
                                      Ipopt::Number *x_u, Ipopt::Index m,
                                      Ipopt::Number *g_l, Ipopt::Number *g_u)
{

    // TODO Check if using information on dimensions could be useful here
    double margin = 2*std::max(object_prop[0],std::max(object_prop[1], object_prop[2]));

    // center
    x_l[0]=centroid[0]-margin; x_u[0]=centroid[0]+margin;
    x_l[1]=centroid[1]-margin; x_u[1]=centroid[1]+margin;
    x_l[2]=centroid[2]-margin; x_u[2]=centroid[2]+margin;

    // three angles in this case
    x_l[3]=0; x_u[3]=0;
    x_l[4]=0; x_u[4]=0;
    x_l[5]=-numeric_limits<double>::infinity(); x_u[5]=numeric_limits<double>::infinity();

    return true;
}

/****************************************************************/
bool SuperQuadricVerticalNLP::get_starting_point(Ipopt::Index n, bool init_x,
                                         Ipopt::Number *x, bool init_z,
                                         Ipopt::Number *z_L, Ipopt::Number *z_U,
                                         Ipopt::Index m, bool init_lambda,
                                         Ipopt::Number *lambda)
{
    x[0]=centroid[0];
    x[1]=centroid[1];
    x[2]=centroid[2];

    x[3]=0;
    x[4]=0;
    x[5]=initial_angle;

    return true;
}

/****************************************************************/
bool SuperQuadricVerticalNLP::eval_f(Ipopt::Index n, const Ipopt::Number *x,
                             bool new_x, Ipopt::Number &obj_value)
{
    Vector c(3),s(3);
    c[0]=x[0];
    c[1]=x[1];
    c[2]=x[2];
    const double &y=x[5];
    s[0]=object_prop[0];
    s[1]=object_prop[1];
    s[2]=object_prop[2];
    const double &e1=object_prop[3];
    const double &e2=object_prop[4];

    Vector angles(3,0.0);
    angles[2] = y;

    Matrix T=rpy2dcm(angles);
    T.setSubcol(c,0,3);
    T=SE3inv(T);

    obj_value=0.0;
    Vector p1(4,1.0);

    for (auto &p:points)
    {
        p1.setSubvector(0,p);
        p1=T*p1;
        double tx=pow(abs( p1[0]/s[0]), 2.0/e2);
        double ty=pow(abs( p1[1]/s[1]), 2.0/e2);
        double tz=pow(abs( p1[2]/s[2]), 2.0/e1);
        double F1=pow(pow( tx+ty, e2/e1) + tz, e1)-1.0;
        obj_value += F1 * F1;
    }

    obj_value*=(s[0]*s[1]*s[2])/points.size();

    return true;
}

/****************************************************************/
bool SuperQuadricVerticalNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number *x,
                                  bool new_x, Ipopt::Number *grad_f)
{
    Vector s(3);
    s[0]=object_prop[0];
    s[1]=object_prop[1];
    s[2]=object_prop[2];
    const double &e1=object_prop[3];
    const double &e2=object_prop[4];

    for (Ipopt::Index i=0; i<n; i++)
        grad_f[i]=0.0;

    if (analytic)
    {
        Vector c(3);
        c[0]=x[0];
        c[1]=x[1];
        c[2]=x[2];
        const double &y=x[5];

        double cy = cos(y);
        double sy = sin(y);

        Matrix Rz(3,3);
        Rz[2][2] = 1;
        Rz[0][0] = cy;
        Rz[1][1] = cy;
        Rz[0][1] = -sy;
        Rz[1][0] = sy;

        Matrix dRz(3,3);
        dRz[0][0] = -sy;
        dRz[1][1] = -sy;
        dRz[0][1] = -cy;
        dRz[1][0] = cy;

        Matrix invRz = Rz.transposed();

        Matrix dX_dT = -1.0 * invRz;
        Matrix dX_dO_y = dRz.transposed();

        Vector invT = -1.0 * invRz * c;

        for (auto &p:points)
        {
            Vector X = invRz*p+invT;

            double tx = pow( fabs( X[0]/s[0] ), 2.0/e2);
            double ty = pow( fabs( X[1]/s[1] ), 2.0/e2);
            double tz = pow( fabs( X[2]/s[2] ), 2.0/e1);

            double tmp_1 = pow(tx+ty, e2/e1) + tz;
            double F1 = pow(tmp_1, e1) - 1.0;

            double tmp_2 = pow(tmp_1, e1-1.0);
            double tmp_3 = pow(tx+ty, e2/e1-1.0);

            double dF_dX0 = 2 * tmp_2 * tmp_3 * sign(X[0]) * pow(fabs(X[0]/s[0]), 2.0/e2-1.0) / s[0];
            double dF_dX1 = 2 * tmp_2 * tmp_3 * sign(X[1]) * pow(fabs(X[1]/s[1]), 2.0/e2-1.0) / s[1];
            double dF_dX2 = 2 * tmp_2 * sign(X[2]) * pow(fabs(X[2]/s[2]), 2.0/e1-1.0) / s[2];

            grad_f[0] += F1 * (dF_dX0 * dX_dT[0][0] + dF_dX1 * dX_dT[1][0] + dF_dX2 * dX_dT[2][0]);
            grad_f[1] += F1 * (dF_dX0 * dX_dT[0][1] + dF_dX1 * dX_dT[1][1] + dF_dX2 * dX_dT[2][1]);
            grad_f[2] += F1 * (dF_dX0 * dX_dT[0][2] + dF_dX1 * dX_dT[1][2] + dF_dX2 * dX_dT[2][2]);

            double p0 = p[0] - c[0];
            double p1 = p[1] - c[1];
            double p2 = p[2] - c[2];

            double dX_dO[3];
            dX_dO[0] = dX_dO_y[0][0] * p0 + dX_dO_y[0][1] * p1 + dX_dO_y[0][2] * p2;
            dX_dO[1] = dX_dO_y[1][0] * p0 + dX_dO_y[1][1] * p1 + dX_dO_y[1][2] * p2;
            dX_dO[2] = dX_dO_y[2][0] * p0 + dX_dO_y[2][1] * p1 + dX_dO_y[2][2] * p2;

            grad_f[5] += F1 * (dF_dX0 * dX_dO[0] + dF_dX1 * dX_dO[1] + dF_dX2 * dX_dO[2]);
        }

        double coeff = 2.0 * s[0] * s[1] * s[2] / points.size();
        for (Ipopt::Index i=0; i<n; i++)
            grad_f[i] *= coeff;

        return true;
    }
    else
    {
         Vector x_tmp(n, 0.0);
         double grad_p, grad_n;
         double eps = 1e-8;

         for (Ipopt::Index j = 0; j < n; j++)
             x_tmp[j] = x[j];

         for (Ipopt::Index j = 0; j < n; j++)
         {
             x_tmp[j] = x[j] + eps;

             Vector c(3);
             c[0]=x_tmp[0];
             c[1]=x_tmp[1];
             c[2]=x_tmp[2];
             double y=x_tmp[5];

             Vector angles(3,0.0);
             angles[2] = y;

             Matrix T=rpy2dcm(angles);
             T.setSubcol(c,0,3);
             T=SE3inv(T);

             double obj_value=0.0;
             Vector p1(4,1.0);
             for (auto &p:points)
             {
                 p1.setSubvector(0,p);
                 p1=T*p1;
                 double tx=pow(abs( p1[0]/s[0]), 2.0/e2);
                 double ty=pow(abs( p1[1]/s[1]), 2.0/e2);
                 double tz=pow(abs( p1[2]/s[2]), 2.0/e1);
                 double F1=pow(pow( tx+ty, e2/e1) + tz, e1)-1.0;
                 obj_value += F1 * F1;
             }

             grad_p = s[0]*s[1]*s[2] * obj_value/points.size();

             x_tmp[j] = x[j] - eps;

             c[0]=x_tmp[0];
             c[1]=x_tmp[1];
             c[2]=x_tmp[2];
             y=x_tmp[5];

             angles[2] = y;

             T=rpy2dcm(angles);
             T.setSubcol(c,0,3);
             T=SE3inv(T);

             obj_value=0.0;
             for (auto &p:points)
             {
                 p1.setSubvector(0,p);
                 p1=T*p1;
                 double tx=pow(abs( p1[0]/s[0]), 2.0/e2);
                 double ty=pow(abs( p1[1]/s[1]), 2.0/e2);
                 double tz=pow(abs( p1[2]/s[2]), 2.0/e1);
                 double F1=pow(pow( tx+ty, e2/e1) + tz, e1)-1.0;
                 obj_value += F1 * F1;
             }
             grad_n = s[0]*s[1]*s[2] * obj_value/points.size();

             grad_f[j] = (grad_p-grad_n)/(2*eps);
         }

         return true;
    }

}

/****************************************************************/
bool SuperQuadricVerticalNLP::eval_g(Ipopt::Index n, const Ipopt::Number *x,
                             bool new_x, Ipopt::Index m, Ipopt::Number *g)
{
    return true;
}

/****************************************************************/
bool SuperQuadricVerticalNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x,
                                 bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                                 Ipopt::Index *iRow, Ipopt::Index *jCol,
                                 Ipopt::Number *values)
{
    return true;
}

/****************************************************************/
bool SuperQuadricVerticalNLP::eval_h(Ipopt::Index n, const Ipopt::Number *x,
                             bool new_x, Ipopt::Number obj_factor,
                             Ipopt::Index m, const Ipopt::Number *lambda,
                             bool new_lambda, Ipopt::Index nele_hess,
                             Ipopt::Index *iRow, Ipopt::Index *jCol,
                             Ipopt::Number *values)
{
    return true;
}

/****************************************************************/
void SuperQuadricVerticalNLP::finalize_solution(Ipopt::SolverReturn status,
                                        Ipopt::Index n, const Ipopt::Number *x,
                                        const Ipopt::Number *z_L,
                                        const Ipopt::Number *z_U,
                                        Ipopt::Index m, const Ipopt::Number *g,
                                        const Ipopt::Number *lambda,
                                        Ipopt::Number obj_value,
                                        const Ipopt::IpoptData *ip_data,
                                        Ipopt::IpoptCalculatedQuantities *ip_cq)
{
    int superq_n = 11;
    result.resize(superq_n);
    for (Ipopt::Index i=0; i<3; i++)
        result[i]=x[i];
    for (Ipopt::Index i=3; i<5; i++)
        result[i]=0.0;
    result[5]=x[5]*180.0/M_PI;
    for (Ipopt::Index i=0; i<3; i++)
        result[i+6]=object_prop[i];
    for (Ipopt::Index i=3; i<5; i++)
        result[i+6]=object_prop[i];
}

/****************************************************************/
SuperQuadricVerticalNLP::SuperQuadricVerticalNLP(const vector<Vector> &points_,
                                 const Vector &object_prop_,
                                 bool analytic_) :
                                 points(points_), object_prop(object_prop_),
                                 analytic(analytic_)
{
    bounds.resize(3,2);
    bounds(0,0)=bounds(1,0)=bounds(2,0)=numeric_limits<double>::infinity();
    bounds(0,1)=bounds(1,1)=bounds(2,1)=-numeric_limits<double>::infinity();

    for (auto &p:points)
    {
        if (p[0]<bounds(0,0))
            bounds(0,0)=p[0];
        if (p[0]>bounds(0,1))
            bounds(0,1)=p[0];

        if (p[1]<bounds(1,0))
            bounds(1,0)=p[1];
        if (p[1]>bounds(1,1))
            bounds(1,1)=p[1];

        if (p[2]<bounds(2,0))
            bounds(2,0)=p[2];
        if (p[2]>bounds(2,1))
            bounds(2,1)=p[2];
    }

    // Compute centroid of point cloud
    centroid.resize(3,0.0);
    for (unsigned int i=0; i<centroid.length(); i++)
        centroid[i]=0.5*(bounds(i,0)+bounds(i,1));

    // Compute orientation of point cloud
    Vector mean(3,0.0);
    for (auto& point: points)
        mean += point;
    mean *= 1.0/points.size();

    Matrix M=zeros(3,3);
    Matrix R(3,3);
    Matrix u(3,3);
    Matrix v(3,3);

    Vector s(3,0.0);
    Vector n(3,0.0);
    Vector o(3,0.0);
    Vector a(3,0.0);

    for (auto& point: points)
    {
        double x = point(0)-mean(0);
        double y = point(1)-mean(1);
        double z = point(2)-mean(2);
        M(0,0) += x*x;
        M(0,1) += x*y;
        M(0,2) += x*z;
        M(1,1) += y*y;
        M(2,2) += z*z;
        M(1,2) += y*z;
    }

    M(0,0) = M(0,0)/points.size();
    M(0,1) = M(0,1)/points.size();
    M(0,2) = M(0,2)/points.size();
    M(1,1) = M(1,1)/points.size();
    M(2,2) = M(2,2)/points.size();
    M(1,2) = M(1,2)/points.size();

    M(1,0) = M(0,1);
    M(2,0) = M(0,2);
    M(2,1) = M(1,2);

    SVDJacobi(M,u,s,v);
    n=u.getCol(0);
    o=u.getCol(1);
    a=u.getCol(2);

    if(object_prop[0]>object_prop[1] && object_prop[0]>object_prop[2])
    {
        R.setCol(0,n);
        if(object_prop[1]>object_prop[2])
        {
            R.setCol(1,o);
            R.setCol(2,a);
        }
        else
        {
            R.setCol(2,o);
            R.setCol(1,a);
        }
    }
    else if(object_prop[1]>object_prop[2])
    {
        R.setCol(1,n);
        if(object_prop[0]>object_prop[2])
        {
            R.setCol(0,o);
            R.setCol(2,a);
        }
        else
        {
            R.setCol(2,o);
            R.setCol(0,a);
        }
    }
    else
    {
        R.setCol(2,n);
        if(object_prop[0]>object_prop[1])
        {
            R.setCol(0,o);
            R.setCol(1,a);
        }
        else
        {
            R.setCol(1,o);
            R.setCol(0,a);
        }
    }

    if(det(R)<0)
        R.setCol(0, -1.0*R.getCol(0));

    // Relign with vertical

    Vector z_pc = R.getCol(2);
    Vector z(3,0.0);
    z[2] = 1.0;
    Vector tu = cross(z_pc, z);
    double sin_t = norm(tu);
    if (sin_t>std::numeric_limits<double>::epsilon())
    {
        double cos_t = dot(z_pc, z);
        double theta = atan2(sin_t, cos_t);
        Vector o(4);
        o.setSubvector(0, 1.0/sin_t * tu);
        o[3] = theta;
        R = axis2dcm(o).submatrix(0,2, 0,2) * R;
    }

    initial_angle = atan2(R[2][2]*0.5*(R[1][0]-R[0][1]), R[2][2]*0.5*(R[0][0]+R[1][1]));
}

/****************************************************************/
Vector SuperQuadricVerticalNLP::get_result() const
{
    return result;
}
