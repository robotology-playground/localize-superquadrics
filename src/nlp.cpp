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
    // TEST
    //m=1;
    //nnz_jac_g=2;
    //nnz_h_lag=0;
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
    Vector margin(3);
    margin[0]=0.25*(bounds(0,1)-bounds(0,0));
    margin[1]=0.25*(bounds(1,1)-bounds(1,0));
    margin[2]=0.25*(bounds(2,1)-bounds(2,0));

    // center
    x_l[0]=bounds(0,0)+margin[0]; x_u[0]=bounds(0,1)-margin[0];
    x_l[1]=bounds(1,0)+margin[1]; x_u[1]=bounds(1,1)-margin[1];
    x_l[2]=bounds(2,0)+margin[2]; x_u[2]=bounds(2,1)-margin[2];
    // three angles in this case
    x_l[3]=0.0; x_u[3]=2.0*M_PI;
    x_l[4]=0.0; x_u[4]=M_PI;
    x_l[5]=0.0; x_u[5]=2.0*M_PI;
    // limit on z-min
    // TEST
    //g_l[0]=bounds(2,0); g_u[0]=numeric_limits<double>::infinity();
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
    x[3]=0.0;
    x[4]=0.0;
    x[5]=0.0;
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
        //double penalty=(F1<0.0?inside_penalty:1.0);   // No penalty in this case
        obj_value += F1 * F1; //*penalty;
    }

    obj_value*=(s[0]*s[1]*s[2])/points.size();

    return true;
}

/****************************************************************/
bool SuperQuadricNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number *x,
                                  bool new_x, Ipopt::Number *grad_f)
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

    for (Ipopt::Index i=0; i<n; i++)
        grad_f[i]=0.0;

    double coeff = s[0] * s[1] * s[2];
    Vector p1(4,1.0);

    for (auto &point:points)
    {
        p1.setSubvector(0,point);
        p1 = T * p1;

        double tx = pow( abs(p1[0]/s[0]), 2.0/e2);
        double ty = pow( abs(p1[1]/s[1]), 2.0/e2);
        double tz = pow( abs(p1[2]/s[2]), 2.0/e1);
        double F1 = pow( pow( tx + ty, e2/e1) + tz, e1)-1.0;

        double tmp1 = 2.0 * coeff * F1;

        double t19 = c[0] * cos(p) * cos(y);
        double t18 = sin(r) * sin(y) + cos(r) * cos(y) * sin(p);
        double t17 = cos(r) * sin(y) - cos(y) * sin(p) * sin(r);
        double t16 = cos(r) * cos(y);
        double t15 = cos(y) * sin(r);
        double t14 = c[0] * cos(p) * sin(y);
        double t13 = p1[2] * cos(p) * sin(r);
        double t12 = cos(r) * sin(p) * sin(y);
        double t11 = sin(p) * sin(r) * sin(y);
        double t10 = p1[2] * cos(p) * cos(r);
        double t9 = c[1] * t17 - c[2] * t18 - p1[2] * sin(p) - t19 + p1[0] * cos(p) * cos(y) + p1[1] * cos(p) * sin(y);
        double t8 = c[1] * (t16 + t11) - c[2] * (t15 - t12) + p1[0] * t17 - p1[1] * (t16 + t11) + t14 - t13;
        double t7 = p1[0] * t18 - p1[1] * (t15 - t12) + c[0] * sin(p) - c[2] * cos(p) * cos(r) - c[1] * cos(p) * sin(r) + t10;

        double t6 = sign(t8) * pow( abs(t8)/s[1], 2.0/e2) + sign(t9) * pow( abs(t9)/s[0], 2.0/e2);
        double t5 = pow( abs(t7)/s[2], 2.0/e1 - 1);
        double t4 = pow( abs(t9)/s[0], 2.0/e2 - 1);
        double t3 = pow( abs(t8)/s[1], 2.0/e2 - 1);
        double t2 = pow( t6, e2/e1 - 1);
        double t1 = pow( ( sign(t6) * pow( t6, e2/e1) + sign(t7) * pow( abs(t7)/s[2], 2.0/e1) ), e1 -1);

        double grad0_0 = (2 * sign(t9) * cos(p) * cos(y) * t4) / (e2 * s[0]);
        double grad0_1 = (2 * sign(t8) * cos(p) * sin(y) * t3) / (e2 * s[1]);
        double grad0_2 = (2 * sign(t7) * sin(p) * t5 ) / (e1 * s[2]);

        grad_f[0] += - tmp1 * e1 * t1 * ( e2/e1 * t2 * (grad0_0 - grad0_1) - grad0_2);

        double grad1_0 = (2 * sign(t9) * t4 * t17) / (e2 * s[0]);
        double grad1_1 = (2 * sign(t8) * (t16 + t11) * t3) / (e2 * s[1]);
        double grad1_2 = (2 * sign(t7) * cos(p) * sin(r) * t5) / (e1 * s[2]);

        grad_f[1] += tmp1 * e1 * t1 * ( e2 * t2 / e1 * (grad1_0 + grad1_1) - grad1_2);

        double grad2_0 = (2 * sign(t9) * t4 * t18) / (e2 * s[0]);
        double grad2_1 = (2 * sign(t8) * (t15 - t12) * t3) / (e2 * s[1]);
        double grad2_2 = (2 * sign(t7) * cos(p) * cos(r) * t5) / (e1 * s[2]);

        grad_f[2] += - tmp1 * e1 * t1 * ( e2 * t2/e1 * (grad2_0 + grad2_1) + grad2_2);

        double grad3_0 = (2 * sign(t9) * t4 * (c[1] * t18 + c[2] * t17)) / (e2 * s[0]);
        double grad3_1 = (2 * sign(t8) * t3 * (c[1] * (t15 - t12) + c[2] * (t16 + t11) + p1[0] * t18 - p1[1] * (t15 - t12) + t10)) / (e2 * s[1]);
        double grad3_2 = (2 * sign(t7) * t5 * (p1[1] * (t16 + t11) - p1[0] * t17 + c[1] * cos(p) * cos(r) - c[2] * cos(p) * sin(r) + t13)) / (e1 * s[2]);

        grad_f[3] += -tmp1 * e1 * t1 * (e1 * t2 / e1 * (grad3_0 + grad3_1) + grad3_2);

        double grad4_0 = (2 * sign(t8) * t3 * (c[0] * sin(p) * sin(y) - p1[2] * sin(p) * sin(r) - c[2] *  cos(p) * cos(r) * sin(y) - c[1] * cos(p) * sin(r) * sin(y) +
                                                                                p1[0] * cos(p) * cos(y) * sin(r) + p1[1] * cos(p) * sin(r) * sin(y))) / (e2 * s[1]);
        double grad4_1 = (2 * sign(t9) * t4 * (p1[2] * cos(p) - c[0] * cos(y) * sin(p) + p1[0] * cos(y) * sin(p) + p1[1] * sin(p) * sin(y) + c[2] * cos(p) * cos(r)  * cos(y) +
                                                                                c[1] * cos(p) * cos(y) * sin(r))) /(e2 * s[0]);
        double grad4_2 = (2 * sign(t7) * t5 * (c[0] * cos(p) + c[2] * cos(r) * sin(p) + c[1] * sin(p) * sin(r) - p1[2] * cos(r) * sin(p) + p1[0] * cos(p) * cos(r) * cos(y) +
                                                                                p1[1] * cos(p) * cos(r) * sin(y)))/(e1 * s[2]);

        grad_f[4] += -tmp1 * e1 * t1 * ( ( e2 * t2 * grad4_0 + grad4_1 / e1)  - grad4_2);

        double grad5_0 = (2 * sign(t8) * t3 * (c[2] * t18 - c[1] * t17 - t16 + p1[0] * (t16 + t11) + p1[1] * t17 - t11 + t19)) / ( e2 * s[1]);
        double grad5_1 = (2 * sign(t9) * t4 * ( cos(p) * sin(y) + c[1] * (t16 + t11) - c[2] * (t15 - t12) + t14 + p1[1] * cos(p) * cos(y) - p1[0] * cos(p) * sin(y)))/ (e2 * s[0]);
        double grad5_2 = (2 * sign(t7) * t5 * (p1[0] * (t15 - t12) - t15 + p1[1] * t18 + t12)) / (e1 * s[2]);

        grad_f[5] += tmp1 * e1 * t1 * (( e2 * t2 * (grad5_0 + grad5_1) / e1) + grad5_2);
    }

    for (Ipopt::Index i=0; i<n; i++)
        grad_f[i]/=points.size();

    return true;
}

/****************************************************************/
bool SuperQuadricNLP::eval_g(Ipopt::Index n, const Ipopt::Number *x,
                             bool new_x, Ipopt::Index m, Ipopt::Number *g)
{
    // TEST
    // g[0]=x[2]-x[6];
    // return true;
    return false;
}

/****************************************************************/
bool SuperQuadricNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x,
                                 bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                                 Ipopt::Index *iRow, Ipopt::Index *jCol,
                                 Ipopt::Number *values)
{
    // TEST
    // if (values==nullptr)
    // {
    //     iRow[0]=0; jCol[0]=2;
    //     iRow[1]=0; jCol[1]=6;
    // }
    // else
    // {
    //     values[0]=1.0;
    //     values[1]=-1.0;
    // }
    // return true;
    return false;
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
                                 const double inside_penalty_,
                                 const Vector &object_prop_) :
                                 points(points_), inside_penalty(inside_penalty_),
                                 object_prop(object_prop_)
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

    centroid.resize(3,0.0);
    for (unsigned int i=0; i<centroid.length(); i++)
        centroid[i]=0.5*(bounds(i,0)+bounds(i,1));
}

/****************************************************************/
Vector SuperQuadricNLP::get_result() const
{
    return result;
}
