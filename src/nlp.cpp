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

#include "nlp.h"

using namespace std;
using namespace yarp::sig;
using namespace yarp::math;

/****************************************************************/
bool SuperQuadricNLP::get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
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
bool SuperQuadricNLP::get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l,
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
    x_l[3]=-numeric_limits<double>::infinity(); x_u[3]=numeric_limits<double>::infinity();
    x_l[4]=-numeric_limits<double>::infinity(); x_u[4]=numeric_limits<double>::infinity();
    x_l[5]=-numeric_limits<double>::infinity(); x_u[5]=numeric_limits<double>::infinity();

    return true;
}

/****************************************************************/
bool SuperQuadricNLP::get_starting_point(Ipopt::Index n, bool init_x,
                                         Ipopt::Number *x, bool init_z,
                                         Ipopt::Number *z_L, Ipopt::Number *z_U,
                                         Ipopt::Index m, bool init_lambda,
                                         Ipopt::Number *lambda)
{
    x[0]=centroid[0];
    x[1]=centroid[1];
    x[2]=centroid[2];

    // x[3]=0.0;
    // x[4]=0.0;
    // x[5]=0.0;
    x[3]=initial_angles[0];
    x[4]=initial_angles[1];
    x[5]=initial_angles[2];

    return true;
}

/****************************************************************/
bool SuperQuadricNLP::eval_f(Ipopt::Index n, const Ipopt::Number *x,
                             bool new_x, Ipopt::Number &obj_value)
{
    Vector c(3),s(3);
    c[0]=x[0];
    c[1]=x[1];
    c[2]=x[2];
    const double &r=x[3];
    const double &p=x[4];
    const double &y=x[5];
    s[0]=object_prop[0];
    s[1]=object_prop[1];
    s[2]=object_prop[2];
    const double &e1=object_prop[3];
    const double &e2=object_prop[4];

    Vector angles(3,0.0);
    angles[0] = r;
    angles[1] = p;
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
bool SuperQuadricNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number *x,
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
        const double &r=x[3];
        const double &p=x[4];
        const double &y=x[5];

        double cr = cos(r);
        double cp = cos(p);
        double cy = cos(y);
        double sr = sin(r);
        double sp = sin(p);
        double sy = sin(y);

        double t15 = cr * cy;
        double t14 = sp * sr * sy;
        double t13 = cr * sy - cy * sp * sr;
        double t12 = cy * sr;
        double t11 = cr * sp * sy;
        double t10 = sr * sy + cr * cy * sp;

        for (auto &point:points)
        {
            double t9 = -cp * cy * (c[0]-point[0]) - cp * sy * (c[1]-point[1]) + sp * (c[2]-point[2]);
            double t8 = t13 * (c[0]-point[0]) - (t15 + t14) * (c[1]-point[1]) - cp * sr * (c[2]-point[2]);
            double t7 = t10 * (c[0]-point[0]) - (t12 - t11) * (c[1]-point[1]) + cp * cr * (c[2]-point[2]);

            double t6 = pow( fabs(t9/s[0]), 2.0/e2) + pow( fabs(t8/s[1]), 2.0/e2);

            double t5 = pow( fabs(t9/s[0]), 2.0/e2-1);
            double t4 = pow( fabs(t7/s[2]), 2.0/e1-1);
            double t3 = pow( fabs(t8/s[1]), 2.0/e2-1);

            double t2 = pow( t6, e2/e1-1);
            double t1 = pow( pow( fabs(t7/s[2]), 2.0/e1) + pow(t6, e2/e1), e1-1);

            double F1 = pow( pow( fabs(t7/s[2]), 2.0/e1) + pow(t6, e2/e1), e1) - 1.0;

            double grad0_0 = (2 * sign(t8) * t13 * t3) / (e2 * s[1]);
            double grad0_1 = (2 * sign(t9) * cp * cy * t5) / (e2 * s[0]);
            double grad0_2 = (2 * sign(t7) * t10 * t4) / (e1 * s[2]);

            grad_f[0] += F1 * e1 * t1 * ( e2/e1 * t2 * (grad0_0 - grad0_1) - grad0_2);

            double grad1_0 = (2 * sign(t8) * (t15 + t14) * t3) / (e2 * s[1]);
            double grad1_1 = (2 * sign(t9) * cp * sy * t5) / (e2 * s[0]);
            double grad1_2 = (2 * sign(t7) * (t12 - t11) * t4) / (e1 * s[2]);

            grad_f[1] += -F1 * e1 * t1 * ( e2/e1 * t2 * (grad1_0 + grad1_1) + grad1_2);

            double grad2_0 = (2 * sign(t9) * sp * t5) / (e2 * s[0]);
            double grad2_1 = (2 * sign(t8) * cp * sr * t3) / (e2 * s[1]);
            double grad2_2 = (2 * sign(t7) * cp * cr * t4) / (e1 * s[2]);

            grad_f[2] += F1 * e1 * t1 * ( e2/e1 * t2 * (grad2_0 - grad2_1) + grad2_2);

            double grad3_0 = (2 * sign(t7) * t4 * t8) / (e1 * s[2]);
            double grad3_1 = (2 * sign(t8) * t3 * t2 * t7) / (e1 * s[1]);

            grad_f[3] += F1 * e1 * t1 * (grad3_0 - grad3_1);

            double grad4_0 = 2 * sign(t9) * t5 / (e2 * s[0]) * ( cy * sp * (c[0]-point[0]) + sp * sy * (c[1]-point[1]) + cp * (c[2]-point[2]) );
            double grad4_1 = 2 * sign(t8) * t3 / (e2 * s[1]) * ( - cp * cy * sr * (c[0]-point[0]) - cp * sr * sy * (c[1]-point[1]) + sp * sr * (c[2]-point[2]) );
            double grad4_2 = 2 * sign(t7) * t4 / (e1 * s[2]) * ( - cp * cy * cr * (c[0]-point[0]) - cp * cr * sy * (c[1]-point[1]) + sp * cr * (c[2]-point[2]) );

            grad_f[4] += F1 * e1 * t1 * ( e2/e1 * t2 * (grad4_0 + grad4_1) - grad4_2);

            double grad5_0 = 2 * sign(t8) * t3 / (e2 * s[1]) * ( (t15 + t14) * (c[0]-point[0]) + t13 * (c[1]-point[1]));
            double grad5_1 = 2 * sign(t9) * t5 / (e2 * s[0]) * ( cp * sy * (c[0]-point[0]) - cp * cy * (c[1]-point[1]));
            double grad5_2 = 2 * sign(t7) * t4 / (e1 * s[2]) * ( (t12 - t11) * (c[0]-point[0]) + t10 * (c[1]-point[1]));

            grad_f[5] += F1 * e1 * t1 * (( e2/e1 * t2 * (grad5_0 + grad5_1)) + grad5_2);
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
             double r=x_tmp[3];
             double p=x_tmp[4];
             double y=x_tmp[5];

             Vector angles(3,0.0);
             angles[0] = r;
             angles[1] = p;
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
             r=x_tmp[3];
             p=x_tmp[4];
             y=x_tmp[5];

             angles[0] = r;
             angles[1] = p;
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
bool SuperQuadricNLP::eval_g(Ipopt::Index n, const Ipopt::Number *x,
                             bool new_x, Ipopt::Index m, Ipopt::Number *g)
{
    return true;
}

/****************************************************************/
bool SuperQuadricNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x,
                                 bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                                 Ipopt::Index *iRow, Ipopt::Index *jCol,
                                 Ipopt::Number *values)
{
    return true;
}

/****************************************************************/
bool SuperQuadricNLP::eval_h(Ipopt::Index n, const Ipopt::Number *x,
                             bool new_x, Ipopt::Number obj_factor,
                             Ipopt::Index m, const Ipopt::Number *lambda,
                             bool new_lambda, Ipopt::Index nele_hess,
                             Ipopt::Index *iRow, Ipopt::Index *jCol,
                             Ipopt::Number *values)
{
    return true;
}

/****************************************************************/
void SuperQuadricNLP::finalize_solution(Ipopt::SolverReturn status,
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
    for (Ipopt::Index i=3; i<6; i++)
        result[i]=x[i]*180.0/M_PI;
    for (Ipopt::Index i=0; i<3; i++)
        result[i+6]=object_prop[i];
    for (Ipopt::Index i=3; i<5; i++)
        result[i+6]=object_prop[i];
}

/****************************************************************/
SuperQuadricNLP::SuperQuadricNLP(const vector<Vector> &points_,
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
        M(0,0) = M(0,0) + (point(1)-centroid(1))*(point(1)-centroid(1)) + (point(2)-centroid(2))*(point(2)-centroid(2));
        M(0,1) = M(0,1) - (point(1)-centroid(1))*(point(0)-centroid(0));
        M(0,2) = M(0,2) - (point(2)-centroid(2))*(point(0)-centroid(0));
        M(1,1) = M(1,1) + (point(0)-centroid(0))*(point(0)-centroid(0)) + (point(2)-centroid(2))*(point(2)-centroid(2));
        M(2,2) = M(2,2) + (point(1)-centroid(1))*(point(1)-centroid(1)) + (point(0)-centroid(0))*(point(0)-centroid(0));
        M(1,2) = M(1,2) - (point(2)-centroid(2))*(point(1)-centroid(1));
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

    R.setCol(0,n);
    R.setCol(1,o);
    R.setCol(2,a);

    initial_angles.resize(3,0.0);
    initial_angles = dcm2rpy(R);
}

/****************************************************************/
Vector SuperQuadricNLP::get_result() const
{
    return result;
}
