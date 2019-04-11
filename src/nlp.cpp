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
    n=6; m=1;
    nnz_jac_g=2;
    nnz_h_lag=0;
    index_style=TNLP::C_STYLE;
    return true;
}

/****************************************************************/
bool SuperQuadricNLP::get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l,
                                      Ipopt::Number *x_u, Ipopt::Index m,
                                      Ipopt::Number *g_l, Ipopt::Number *g_u)
{
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
    g_l[0]=bounds(2,0); g_u[0]=numeric_limits<double>::infinity();
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
    const double &a1=x[3];
    const double &a2=x[4];
    const double &a3=x[5];
    s[0]=object_prop[0];
    s[1]=object_prop[1];
    s[2]=object_prop[2];
    const double &e1=object_prop[3];
    const double &e2=object_prop[4];

    // TODO To change this becasue we will estimate position and orientation

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
    const double &a1=x[3];
    const double &a2=x[4];
    const double &a3=x[5];
    s[0]=object_prop[0];
    s[1]=object_prop[1];
    s[2]=object_prop[2];
    const double &e1=object_prop[3];
    const double &e2=object_prop[4];

    // TODO To change this becasue we will estimate position and orientation

    return true;
}

/****************************************************************/
bool SuperQuadricNLP::eval_g(Ipopt::Index n, const Ipopt::Number *x,
                             bool new_x, Ipopt::Index m, Ipopt::Number *g)
{
    g[0]=x[2]-x[6];
    return true;
}

/****************************************************************/
bool SuperQuadricNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x,
                                 bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                                 Ipopt::Index *iRow, Ipopt::Index *jCol,
                                 Ipopt::Number *values)
{
    if (values==nullptr)
    {
        iRow[0]=0; jCol[0]=2;
        iRow[1]=0; jCol[1]=6;
    }
    else
    {
        values[0]=1.0;
        values[1]=-1.0;
    }
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
    // TODO Check if the result is stored properly
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
