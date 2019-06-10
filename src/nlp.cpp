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
Vector SuperQuadricNLP::projectPoint(const yarp::sig::Vector &point)
{
    Vector pixel_lambda(3,0.0);
    Vector pixel(2);

    pixel_lambda=K*T*point;
    pixel[0]=pixel_lambda[0]/pixel_lambda[2];
    pixel[1]=pixel_lambda[1]/pixel_lambda[2];

    return pixel;
}

/****************************************************************/
double SuperQuadricNLP::insideMask(Vector &point)
{
     int nbSides=object_contour.size();
     int on_right=0;
     double dist=0;
     double min_dist=std::numeric_limits<double>::max();
     Vector p1(2);
     Vector p2=object_contour[nbSides-1];
     Vector x(2);
     Vector x1(2);
     Vector x2(2);

     for(size_t i=0; i<nbSides; i++)
     {
         p1=p2;
         p2=object_contour[i];

         //cout << "object_contour " << object_contour[i].toString() << endl;
         // cout << "p1 " << p1.toString() << endl;

         x=p2-p1;
         x*=1.0/norm(x);
         x1=point-p1;
         x2=point-p2;

         if( dot(x1,x)<=0 )
             dist=norm(x1);
         else if( dot(x2,x)>=0 )
             dist=norm(x2);
         else
             dist=fabs(x1[1]*x[0]-x1[0]*x[1]);

         if(dist<min_dist)
         {
             min_dist=dist;
             if(min_dist==0)
                 break;
         }

         if( (p1[1]<=point[1] && p2[1]<=point[1]) ||
            (p1[1]>point[1] && p2[1]>point[1]) ||
            (p1[0]<point[0] && p2[0]<point[0]) )
             continue;

         double signed_dist=x1[1]*x[0]-x1[0]*x[1];
         if(x[1]<0)
             signed_dist=-signed_dist;
         if(signed_dist>0)
             on_right++;
     }

     if(on_right%2==0)
     {
         return -min_dist;
     }
     else
     {
         return min_dist;
     }
}

/****************************************************************/
bool SuperQuadricNLP::get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
                                   Ipopt::Index &nnz_jac_g,
                                   Ipopt::Index &nnz_h_lag,
                                   IndexStyleEnum &index_style)
{
    n=6;
    m=1;
    nnz_jac_g=n*m;
    nnz_h_lag=0;
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

    g_l[0]=g_u[0]=0.0;

    return true;
}

/****************************************************************/
bool SuperQuadricNLP::get_starting_point(Ipopt::Index n, bool init_x,
                                         Ipopt::Number *x, bool init_z,
                                         Ipopt::Number *z_L, Ipopt::Number *z_U,
                                         Ipopt::Index m, bool init_lambda,
                                         Ipopt::Number *lambda)
{
    x[0]=centroid[0] + 0.05;
    x[1]=centroid[1];
    x[2]=centroid[2];
    // x[0]=0.0;
    // x[1]=0.0;
    // x[2]=0.0;
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

    Vector p1(4,1.0);
    double distance=0.0;

    for (auto p: points_superq)
    {
        p1.setSubvector(0,p);
        p1=T*p1;
        Vector pixel=projectPoint(p1);

        double d=-insideMask(pixel);
        if (d > 0)
            distance+=d;
    }

    obj_value = distance/points.size();

    return true;
}

/****************************************************************/
double SuperQuadricNLP::F_v(Vector &x)
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

    Vector p1(4,1.0);
    double distance=0.0;

    for (auto p: points_superq)
    {
        p1.setSubvector(0,p);
        p1=T*p1;
        Vector pixel=projectPoint(p1);
        double d=-insideMask(pixel);

        if (d > 0)
        {
            // cout << "p " << p.toString() << endl;
            // cout << "d " << d << endl;
            distance+=d;
        }
    }

    //cout << "distance " << distance << endl;
    return distance/points.size();
}

/****************************************************************/
bool SuperQuadricNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number *x,
                                  bool new_x, Ipopt::Number *grad_f)
{
     Vector x_tmp(n,0.0);
     double grad_p, grad_n;
     double eps = 1e-8;

     for(Ipopt::Index i = 0;i < n; i++)
        x_tmp(i) = x[i];

     for(Ipopt::Index j = 0;j < n; j++)
     {
         x_tmp(j) += eps;
         grad_p = F_v(x_tmp);

         // cout<< "x_tmp " << x_tmp.toString() << endl;
         // cout<< "grad_p " << grad_p << endl;

         x_tmp(j) -= eps;
         grad_n = F_v(x_tmp);

         // cout<< "x_tmp " << x_tmp.toString() << endl;
         // cout<< "grad_n " << grad_n << endl;

         grad_f[j] = (grad_p-grad_n)/eps;
         cout<< "Gradient value " << grad_f[j] << endl;
     }

    return true;
}

/****************************************************************/
bool SuperQuadricNLP::eval_g(Ipopt::Index n, const Ipopt::Number *x,
                             bool new_x, Ipopt::Index m, Ipopt::Number *g)
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

    g[0]=obj_value*(s[0]*s[1]*s[2])/points.size();

    cout << "Constraint value " << g[0] << endl;

    return true;
}

/****************************************************************/
double SuperQuadricNLP::G_v(Vector &x)
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

     //cout<< "Constraint value " << obj_value*(s[0]*s[1]*s[2])/points.size() << endl;

    return obj_value*(s[0]*s[1]*s[2])/points.size();
}

/****************************************************************/
bool SuperQuadricNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x,
                                 bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                                 Ipopt::Index *iRow, Ipopt::Index *jCol,
                                 Ipopt::Number *values)
{
     double grad_p, grad_n;
     double eps = 1e-6;
     Vector x_tmp(6,0.0);

     if(values != NULL)
     {
         for(Ipopt::Index i = 0;i < n; i++)
            x_tmp(i) = x[i];

         int count = 0;
         for(Ipopt::Index i = 0;i < m; i++)
         {
             for(Ipopt::Index j = 0;j < n; j++)
             {
                 x_tmp(j) = x_tmp(j) + eps;
                 grad_p = G_v(x_tmp);

                 x_tmp(j) = x_tmp(j)-eps;
                 grad_n = G_v(x_tmp);

                 values[count] = (grad_p-grad_n)/(eps);
                 count++;
             }
         }
     }
     else
     {
        for (int j = 0; j < m; j++)
        {
            for (int i = 0; i < n; i++)
            {
                jCol[j*(n) + i] = i;
                iRow[j*(n)+ i] = j;
            }
        }
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
                                 const vector<Vector> &points_superq_,
                                 const vector<Vector> &object_contour_,
                                 const Vector &object_prop_,
                                 const Matrix &K_,
                                 const Matrix &T_,
                                 bool analytic_) :
                                 points(points_), points_superq(points_superq_),
                                 object_contour(object_contour_),
                                 object_prop(object_prop_), K(K_), T(T_),
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
