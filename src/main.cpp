/******************************************************************************
 *                                                                            *
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file main.cpp
 * @authors: Ugo Pattacini <ugo.pattacini@iit.it>
 */

#include <cstdlib>
#include <memory>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>

#include <yarp/os/all.h>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>
#include <yarp/math/Rand.h>

#include <iCub/ctrl/clustering.h>

#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkSuperquadric.h>
#include <vtkUnsignedCharArray.h>
#include <vtkTransform.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleSwitch.h>

#include "nlp.h"
#include "test.h"
#include "src/Localizer_IDL.h"

using namespace std;
using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;
using namespace iCub::ctrl;

Mutex mutex;

/****************************************************************/
class UpdateCommand : public vtkCommand
{
    const bool *closing;

public:
    /****************************************************************/
    vtkTypeMacro(UpdateCommand, vtkCommand);

    /****************************************************************/
    static UpdateCommand *New()
    {
        return new UpdateCommand;
    }

    /****************************************************************/
    UpdateCommand() : closing(nullptr) { }

    /****************************************************************/
    void set_closing(const bool &closing)
    {
        this->closing=&closing;
    }

    /****************************************************************/
    void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId),
                 void *vtkNotUsed(callData))
    {
        LockGuard lg(mutex);
        vtkRenderWindowInteractor* iren=static_cast<vtkRenderWindowInteractor*>(caller);
        if (closing!=nullptr)
        {
            if (*closing)
            {
                iren->GetRenderWindow()->Finalize();
                iren->TerminateApp();
                return;
            }
        }

        iren->GetRenderWindow()->SetWindowName("Find Superquadric");
        iren->Render();
    }
};


/****************************************************************/
class Object
{
protected:
    vtkSmartPointer<vtkPolyDataMapper> vtk_mapper;
    vtkSmartPointer<vtkActor> vtk_actor;

public:
    /****************************************************************/
    vtkSmartPointer<vtkActor> &get_actor()
    {
        return vtk_actor;
    }
};


/****************************************************************/
class Points : public Object
{
protected:
    vtkSmartPointer<vtkPoints> vtk_points;
    vtkSmartPointer<vtkUnsignedCharArray> vtk_colors;
    vtkSmartPointer<vtkPolyData> vtk_polydata;
    vtkSmartPointer<vtkVertexGlyphFilter> vtk_glyphFilter;

public:
    /****************************************************************/
    Points(const vector<Vector> &points, const int point_size)
    {
        vtk_points=vtkSmartPointer<vtkPoints>::New();
        for (size_t i=0; i<points.size(); i++)
            vtk_points->InsertNextPoint(points[i][0],points[i][1],points[i][2]);

        vtk_polydata=vtkSmartPointer<vtkPolyData>::New();
        vtk_polydata->SetPoints(vtk_points);

        vtk_glyphFilter=vtkSmartPointer<vtkVertexGlyphFilter>::New();
        vtk_glyphFilter->SetInputData(vtk_polydata);
        vtk_glyphFilter->Update();

        vtk_mapper=vtkSmartPointer<vtkPolyDataMapper>::New();
        vtk_mapper->SetInputConnection(vtk_glyphFilter->GetOutputPort());

        vtk_actor=vtkSmartPointer<vtkActor>::New();
        vtk_actor->SetMapper(vtk_mapper);
        vtk_actor->GetProperty()->SetPointSize(point_size);
    }

    /****************************************************************/
    void set_points(const vector<Vector> &points)
    {
        vtk_points=vtkSmartPointer<vtkPoints>::New();
        for (size_t i=0; i<points.size(); i++)
            vtk_points->InsertNextPoint(points[i][0],points[i][1],points[i][2]);

        vtk_polydata->SetPoints(vtk_points);
    }

    /****************************************************************/
    bool set_colors(const vector<vector<unsigned char>> &colors)
    {
        if (colors.size()==vtk_points->GetNumberOfPoints())
        {
            vtk_colors=vtkSmartPointer<vtkUnsignedCharArray>::New();
            vtk_colors->SetNumberOfComponents(3);
            for (size_t i=0; i<colors.size(); i++)
                vtk_colors->InsertNextTypedTuple(colors[i].data());

            vtk_polydata->GetPointData()->SetScalars(vtk_colors);
            return true;
        }
        else
            return false;
    }

    /****************************************************************/
    vtkSmartPointer<vtkPolyData> &get_polydata()
    {
        return vtk_polydata;
    }
};


/****************************************************************/
class Superquadric : public Object
{
protected:
    vtkSmartPointer<vtkSuperquadric> vtk_superquadric;
    vtkSmartPointer<vtkSampleFunction> vtk_sample;
    vtkSmartPointer<vtkContourFilter> vtk_contours;
    vtkSmartPointer<vtkTransform> vtk_transform;

public:
    /****************************************************************/
    Superquadric(const Vector &r, const vector<double> &color,
                 const double opacity)
    {
        double bx=r[6];
        double by=r[8];
        double bz=r[7];

        vtk_superquadric=vtkSmartPointer<vtkSuperquadric>::New();
        vtk_superquadric->ToroidalOff();
        vtk_superquadric->SetSize(1.0);
        vtk_superquadric->SetCenter(zeros(3).data());

        vtk_superquadric->SetScale(bx,by,bz);
        vtk_superquadric->SetPhiRoundness(r[9]); // roundness along model z axis (vtk y axis)
        vtk_superquadric->SetThetaRoundness(r[10]); // common roundness along model x and y axes (vtk x and

        vtk_sample=vtkSmartPointer<vtkSampleFunction>::New();
        vtk_sample->SetSampleDimensions(50,50,50);
        vtk_sample->SetImplicitFunction(vtk_superquadric);
        vtk_sample->SetModelBounds(-bx,bx,-by,by,-bz,bz);

        // The isosurface is defined at 0.0 as specified in
        // https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkSuperquadric.cxx
        vtk_contours=vtkSmartPointer<vtkContourFilter>::New();
        vtk_contours->SetInputConnection(vtk_sample->GetOutputPort());
        vtk_contours->GenerateValues(1,0.0,0.0);

        vtk_mapper=vtkSmartPointer<vtkPolyDataMapper>::New();
        vtk_mapper->SetInputConnection(vtk_contours->GetOutputPort());
        vtk_mapper->ScalarVisibilityOff();

        vtk_actor=vtkSmartPointer<vtkActor>::New();
        vtk_actor->SetMapper(vtk_mapper);
        vtk_actor->GetProperty()->SetColor(color[0],color[1],color[2]);
        vtk_actor->GetProperty()->SetOpacity(opacity);

        vtk_transform=vtkSmartPointer<vtkTransform>::New();
        vtk_transform->Translate(r.subVector(0,2).data());

        Vector axisangle = dcm2axis(rpy2dcm(r.subVector(3,5)* M_PI/180.0));
        vtk_transform->RotateWXYZ((180.0/M_PI)*axisangle[3], axisangle.subVector(0,2).data());
        vtk_transform->RotateX(-90.0); // rotate to invert y and z
        vtk_actor->SetUserTransform(vtk_transform);
    }

    /****************************************************************/
    void set_parameters(const Vector &r)
    {
        // Note: roundness parameter for axes x and y is shared in SQ model,
        //       but VTK shares axes x and z (ThetaRoundness).
        //       To get a good display, directions of axes y and z need to be swapped
        //       => parameters for y and z are inverted and a rotation of -90 degrees around x is added

        double bx=r[6];
        double by=r[8];
        double bz=r[7];

        vtk_superquadric->SetScale(bx,by,bz);
        vtk_superquadric->SetPhiRoundness(r[9]); // roundness along model z axis (vtk y axis)
        vtk_superquadric->SetThetaRoundness(r[10]); // common roundness along model x and y axes (vtk x and z axes)

        vtk_sample->SetModelBounds(-bx,bx,-by,by,-bz,bz);

        Vector axisangle = dcm2axis(rpy2dcm(r.subVector(3,5) * M_PI/180.0));

        vtk_transform->Identity();
        vtk_transform->Translate(r.subVector(0,2).data());
        vtk_transform->RotateWXYZ((180.0/M_PI)*axisangle[3], axisangle.subVector(0,2).data());
        vtk_transform->RotateX(-90.0); // rotate to invert y and z
        vtk_actor->SetUserTransform(vtk_transform);
    }
};


/****************************************************************/
class Localizer : public RFModule, Localizer_IDL
{
    Bottle outliersRemovalOptions;
    unsigned int uniform_sample;
    double random_sample;
    bool from_file;
    bool test_derivative;
    bool viewer_enabled;
    bool analytic_gradient;
    bool closing;

    // New
    int num_vis;
    string object_name;
    Vector object_prop;
    PointCloudXYZRGBA point_cloud;

    // The rpc port now is thrifted and used also for
    // localizing the superquadric
    RpcServer rpcService;

    // This rpc to communicate with OPC
    RpcClient rpcOPC;

    vector<Vector> all_points,in_points,out_points,dwn_points;
    vector<vector<unsigned char>> all_colors;

    vector<unique_ptr<Points>> vtk_all_points,vtk_out_points,vtk_dwn_points;
    vector<unique_ptr<Superquadric>> vtk_superquadrics;

    vtkSmartPointer<vtkRenderer> vtk_renderer;
    vtkSmartPointer<vtkRenderWindow> vtk_renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> vtk_renderWindowInteractor;
    vtkSmartPointer<vtkAxesActor> vtk_axes;
    vtkSmartPointer<vtkOrientationMarkerWidget> vtk_widget;
    vtkSmartPointer<vtkCamera> vtk_camera;
    vtkSmartPointer<vtkInteractorStyleSwitch> vtk_style;
    vtkSmartPointer<UpdateCommand> vtk_updateCallback;

    /************************************************************************/
    bool attach(RpcServer &source)
    {
      return this->yarp().attachAsServer(source);
    }

    /****************************************************************/
    void removeOutliers()
    {
        double t0=Time::now();
        if (outliersRemovalOptions.size()>=2)
        {
            double radius=outliersRemovalOptions.get(0).asDouble();
            int minpts=outliersRemovalOptions.get(1).asInt();

            Property options;
            options.put("epsilon",radius);
            options.put("minpts",minpts);

            DBSCAN dbscan;
            map<size_t,set<size_t>> clusters=dbscan.cluster(all_points,options);

            size_t largest_class; size_t largest_size=0;
            for (auto it=begin(clusters); it!=end(clusters); it++)
            {
                if (it->second.size()>largest_size)
                {
                    largest_size=it->second.size();
                    largest_class=it->first;
                }
            }

            auto &c=clusters[largest_class];
            for (size_t i=0; i<all_points.size(); i++)
            {
                if (c.find(i)==end(c))
                    out_points.push_back(all_points[i]);
                else
                    in_points.push_back(all_points[i]);
            }
        }
        else
            in_points=all_points;

        double t1=Time::now();
        yInfo()<<out_points.size()<<"outliers removed out of"
               <<all_points.size()<<"points in"<<t1-t0<<"[s]";
    }

    /****************************************************************/
    void sampleInliers()
    {
        double t0=Time::now();
        if (random_sample>=1.0)
        {
            unsigned int cnt=0;
            for (auto &p:in_points)
            {
                if ((cnt++%uniform_sample)==0)
                    dwn_points.push_back(p);
            }
        }
        else
        {
            set<unsigned int> idx;
            while (idx.size()<(size_t)(random_sample*in_points.size()))
            {
                unsigned int i=(unsigned int)(Rand::scalar(0.0,1.0)*in_points.size());
                if (idx.find(i)==idx.end())
                {
                    dwn_points.push_back(in_points[i]);
                    idx.insert(i);
                }
            }
        }

        double t1=Time::now();
        yInfo()<<dwn_points.size()<<"samples out of"
               <<in_points.size()<<"inliers in"<<t1-t0<<"[s]";
    }

    /****************************************************************/
    Vector localizeSuperquadric() const
    {
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app=new Ipopt::IpoptApplication;
        app->Options()->SetNumericValue("tol",1e-6);
        app->Options()->SetIntegerValue("acceptable_iter",0);
        app->Options()->SetStringValue("mu_strategy","adaptive");
        app->Options()->SetStringValue("nlp_scaling_method","gradient-based");
        app->Options()->SetIntegerValue("max_iter", 1000);
        app->Options()->SetStringValue("hessian_approximation","limited-memory");
        app->Options()->SetStringValue("derivative_test",test_derivative?"first-order":"none");
        app->Options()->SetIntegerValue("print_level",test_derivative?5:0);
        app->Initialize();

        double t0=Time::now();
        Ipopt::SmartPtr<SuperQuadricNLP> nlp=new SuperQuadricNLP(dwn_points, object_prop, analytic_gradient);
        Ipopt::ApplicationReturnStatus status=app->OptimizeTNLP(GetRawPtr(nlp));
        double t1=Time::now();

        Vector r=nlp->get_result();
        yInfo()<<"center   = ("<<r.subVector(0,2).toString(3,3)<<")";
        yInfo()<<"orientation    ="<<r.subVector(3,5).toString(3,3)<<"[deg]";
        yInfo()<<"found in ="<<t1-t0<<"[s]";

        return r;
    }

    /****************************************************************/
    bool configure(ResourceFinder &rf) override
    {
        Rand::init();

        object_prop.resize(5);

        num_vis = 5;

        from_file=rf.check("file");
        if (from_file)
        {
            string file=rf.find("file").asString();
            Bottle *list=rf.find("object_prop").asList();
            if(list)
            {
                if(list->size() == 5)
                {
                    for(int i = 0; i < 5; i++) object_prop[i] = list->get(i).asDouble();
                }
                else
                {
                    yError() <<  "Invalid object shape and dim dimension in config. Should be 5.";
                    return false;
                }
            }

            ifstream fin(file.c_str());
            if (!fin.is_open())
            {
                yError()<<"Unable to open file \""<<file<<"\"";
                return false;
            }

            Vector p(3);
            vector<unsigned int> c_(3);
            vector<unsigned char> c(3);

            string line;
            while (getline(fin,line))
            {
                istringstream iss(line);
                if (!(iss>>p[0]>>p[1]>>p[2]))
                    break;
                all_points.push_back(p);

                fill(c_.begin(),c_.end(),120);
                iss>>c_[0]>>c_[1]>>c_[2];
                c[0]=(unsigned char)c_[0];
                c[1]=(unsigned char)c_[1];
                c[2]=(unsigned char)c_[2];
                all_colors.push_back(c);
            }
        }
        else
        {
            rpcService.open("/localize-superquadric/service:rpc");
            attach(rpcService);

            rpcOPC.open("/localize-superquadric/opc:rpc");
        }

        if (rf.check("remove-outliers"))
            if (const Bottle *ptr=rf.find("remove-outliers").asList())
                outliersRemovalOptions=*ptr;

        uniform_sample=(unsigned int)rf.check("uniform-sample",Value(1)).asInt();
        random_sample=rf.check("random-sample",Value(1.0)).asDouble();
        analytic_gradient=rf.check("analytic-gradient");
        test_derivative=rf.check("test-derivative");
        viewer_enabled=!rf.check("disable-viewer");

        vector<double> color={0.0,0.3,0.6};
        if (rf.check("color"))
        {
            if (const Bottle *ptr=rf.find("color").asList())
            {
                size_t len=std::min(color.size(),ptr->size());
                for (size_t i=0; i<len; i++)
                    color[i]=ptr->get(i).asDouble();
            }
        }

        double opacity=rf.check("opacity",Value(0.25)).asDouble();

        vector<double> backgroundColor={0.7,0.7,0.7};
        if (rf.check("background-color"))
        {
            if (const Bottle *ptr=rf.find("background-color").asList())
            {
                size_t len=std::min(backgroundColor.size(),ptr->size());
                for (size_t i=0; i<len; i++)
                    backgroundColor[i]=ptr->get(i).asDouble();
            }
        }

        removeOutliers();
        sampleInliers();

        Vector r(11,0.0);
        if (dwn_points.size()>0)
            r=localizeSuperquadric();

        vtk_renderer=vtkSmartPointer<vtkRenderer>::New();
        vtk_renderWindow=vtkSmartPointer<vtkRenderWindow>::New();
        vtk_renderWindow->SetSize(600,600);
        vtk_renderWindow->AddRenderer(vtk_renderer);
        vtk_renderWindowInteractor=vtkSmartPointer<vtkRenderWindowInteractor>::New();
        vtk_renderWindowInteractor->SetRenderWindow(vtk_renderWindow);

        for (int i = 0; i < num_vis; i++)
        {
            //Vector r(11,0.0);
            vtk_superquadrics.push_back(unique_ptr<Superquadric>(new Superquadric(r, color,opacity)));
            vtk_renderer->AddActor(vtk_superquadrics[i]->get_actor());

            vtk_all_points.push_back(unique_ptr<Points>(new Points(all_points,2)));
            vtk_out_points.push_back(unique_ptr<Points>(new Points(out_points,4)));
            vtk_dwn_points.push_back(unique_ptr<Points>(new Points(dwn_points,1)));

            vtk_all_points[i]->set_colors(all_colors);
            vtk_out_points[i]->get_actor()->GetProperty()->SetColor(1.0,0.0,0.0);
            vtk_dwn_points[i]->get_actor()->GetProperty()->SetColor(1.0,1.0,0.0);

            vtk_renderer->AddActor(vtk_all_points[i]->get_actor());
            vtk_renderer->AddActor(vtk_out_points[i]->get_actor());
            if (dwn_points.size()!=in_points.size())
                vtk_renderer->AddActor(vtk_dwn_points[i]->get_actor());

        }

        vtk_renderer->SetBackground(backgroundColor.data());

        vtk_axes=vtkSmartPointer<vtkAxesActor>::New();
        vtk_widget=vtkSmartPointer<vtkOrientationMarkerWidget>::New();
        vtk_widget->SetOutlineColor(0.9300,0.5700,0.1300);
        vtk_widget->SetOrientationMarker(vtk_axes);
        vtk_widget->SetInteractor(vtk_renderWindowInteractor);
        vtk_widget->SetViewport(0.0,0.0,0.2,0.2);
        vtk_widget->SetEnabled(1);
        vtk_widget->InteractiveOn();

        vector<double> bounds(6),centroid(3);
        vtk_all_points[0]->get_polydata()->GetBounds(bounds.data());
        for (size_t i=0; i<centroid.size(); i++)
            centroid[i]=0.5*(bounds[i<<1]+bounds[(i<<1)+1]);

        vtk_camera=vtkSmartPointer<vtkCamera>::New();
        vtk_camera->SetPosition(centroid[0]+1.0,centroid[1],centroid[2]+0.5);
        vtk_camera->SetFocalPoint(centroid.data());
        vtk_camera->SetViewUp(0.0,0.0,1.0);
        vtk_renderer->SetActiveCamera(vtk_camera);

        vtk_style=vtkSmartPointer<vtkInteractorStyleSwitch>::New();
        vtk_style->SetCurrentStyleToTrackballCamera();
        vtk_renderWindowInteractor->SetInteractorStyle(vtk_style);

        if (viewer_enabled)
        {
            vtk_renderWindowInteractor->Initialize();
            vtk_renderWindowInteractor->CreateRepeatingTimer(10);

            vtk_updateCallback=vtkSmartPointer<UpdateCommand>::New();
            vtk_updateCallback->set_closing(closing);
            vtk_renderWindowInteractor->AddObserver(vtkCommand::TimerEvent,vtk_updateCallback);
            vtk_renderWindowInteractor->Start();
        }

        return true;
    }

    /****************************************************************/
    double getPeriod() override
    {
        return 1.0;
    }

    /****************************************************************/
    bool updateModule() override
    {
        return (!from_file && !viewer_enabled);
    }

    /****************************************************************/
    bool reset_visualization()
    {
        Vector r(11,0.0);
        all_points.clear();
        all_colors.clear();
        out_points.clear();
        dwn_points.clear();

        for (size_t i = 0 ; i < num_vis ; i++)
        {
            vtk_all_points[i]->set_points(all_points);
            vtk_all_points[i]->set_colors(all_colors);
            vtk_out_points[i]->set_points(out_points);
            vtk_dwn_points[i]->set_points(dwn_points);
            vtk_superquadrics[i]->set_parameters(r);
        }

        return true;
    }

    /****************************************************************/
    Vector localize_superq(const string &object_name, const Bottle &points_bottle, int object_num)
    {
        bool success = point_cloud.fromBottle(points_bottle);
        if (success)
        {
            if (point_cloud.size()>0)
            {
                LockGuard lg(mutex);

                all_points.clear();
                all_colors.clear();
                in_points.clear();
                out_points.clear();
                dwn_points.clear();

                Vector p(3);
                vector<unsigned char> c(3);
                for (int i=0; i<point_cloud.size(); i++)
                {
                    p[0]=point_cloud(i).x;
                    p[1]=point_cloud(i).y;
                    p[2]=point_cloud(i).z;
                    c[0]=point_cloud(i).r;
                    c[1]=point_cloud(i).g;
                    c[2]=point_cloud(i).b;
                    all_points.push_back(p);
                    all_colors.push_back(c);
                }

                // Ask object dimensions and shape to OPC
                Bottle cmd,reply;
                cmd.addVocab(Vocab::encode("ask"));
                Bottle &content=cmd.addList();
                Bottle &cond_1=content.addList();
                cond_1.addString("entity");
                cond_1.addString("==");
                cond_1.addString("object");
                content.addString("&&");
                Bottle &cond_2=content.addList();
                cond_2.addString("name");
                cond_2.addString("==");
                cond_2.addString(object_name);

                rpcOPC.write(cmd, reply);

                if(reply.size()>1)
                {
                    if(reply.get(0).asVocab()==Vocab::encode("ack"))
                    {
                        if (Bottle *b=reply.get(1).asList())
                        {
                            if (Bottle *b1=b->get(1).asList())
                            {
                                cmd.clear();
                                int id=b1->get(0).asInt();
                                cmd.addVocab(Vocab::encode("get"));
                                Bottle &info=cmd.addList();
                                Bottle &info2=info.addList();
                                info2.addString("id");
                                info2.addInt(id);
                                Bottle &info3=info.addList();
                                info3.addString("propSet");
                                Bottle &info4=info3.addList();
                                info4.addString("object_dim_and_shape");
                            }
                            else
                            {
                                yError("no object id provided by OPC!");
                            }
                        }
                        else
                        {
                            yError("uncorrect reply from OPC!");
                        }

                        Bottle reply;
                        if (rpcOPC.write(cmd,reply))
                        {
                            if (reply.size()>1)
                            {
                                if (reply.get(0).asVocab()==Vocab::encode("ack"))
                                {
                                    if (Bottle *b=reply.get(1).asList())
                                    {
                                        if (Bottle *b1=b->find("object_dim_and_shape").asList())
                                        {
                                            object_prop[0]=b1->get(0).asDouble();
                                            object_prop[1]=b1->get(1).asDouble();
                                            object_prop[2]=b1->get(2).asDouble();
                                            object_prop[3]=b1->get(3).asDouble();
                                            object_prop[4]=b1->get(4).asDouble();

                                            yInfo() << "Received object dimension and shape " << object_prop.toString();
                                        }
                                        else
                                        {
                                            yError("object_dim_and_shape field not found in the OPC reply!");
                                        }
                                    }
                                    else
                                    {
                                        yError("uncorrect reply structure received!");
                                    }
                                }
                                else
                                {
                                    yError("Failure in reply for object_dim_and_shape!");
                                }
                            }
                            else
                            {
                                yError("reply size for object_dim_and_shape less than 1!");
                            }
                        }
                        else
                            yError("no reply from second OPC query!");
                    }
                    else
                    {
                        yError("Failure in reply for object id!");
                    }
                }
                else
                {
                    yError("reply size for object id less than 1!");
                }

                removeOutliers();
                sampleInliers();

                Vector r=localizeSuperquadric();

                vtk_all_points[object_num]->set_points(all_points);
                vtk_all_points[object_num]->set_colors(all_colors);
                vtk_out_points[object_num]->set_points(out_points);
                vtk_dwn_points[object_num]->set_points(dwn_points);
                vtk_superquadrics[object_num]->set_parameters(r);

                vtk_camera->SetPosition(r[0]+0.5,r[1],r[2]+0.4);
                vtk_camera->SetFocalPoint(r.subVector(0,2).data());
                vtk_camera->SetViewUp(0.0,0.0,1.0);

                return r;
            }
            else
            {
                yError() << "Point cloud size 0!";

                Vector r(11,0.0);
                return r;
            }

        }
        else
        {
            yError() << "No point cloud bottle received!";

            Vector r(11,0.0);
            return r;
        }
    }

    /****************************************************************/
    bool set_remove_outliers(const Bottle &params)
    {
        if (params.size()==2)
        {
            outliersRemovalOptions.clear();
            double radius=params.get(0).asDouble();
            int minpts=params.get(1).asInt();
            outliersRemovalOptions.addDouble(radius);
            outliersRemovalOptions.addInt(minpts);

            yInfo() << "Set outlier removal options: radius" << outliersRemovalOptions.get(0).asDouble()
                     << "and minpts" <<  outliersRemovalOptions.get(1).asInt();

            return true;
        }
        else
            return false;
    }

    /****************************************************************/
    bool set_uniform_sample(const string &input)
    {
        uniform_sample=(input=="on");

        yInfo() << " New value for uniform_sample " << uniform_sample;

        return true;
    }

    /****************************************************************/
    bool set_random_sample(const string &input)
    {
        random_sample=(input=="on");

        yInfo() << " New value for random_sample " << random_sample;

        return true;
    }

    /****************************************************************/
    bool set_analytic_gradient(const string &input)
    {
        analytic_gradient=(input=="on");

        yInfo() << " New value for analytic_gradient " << analytic_gradient;

        return true;
    }

    /****************************************************************/
    bool interruptModule() override
    {
        closing=true;
        return true;
    }

    /****************************************************************/
    bool close() override
    {
        if (rpcOPC.asPort().isOpen())
            rpcOPC.close();
        if (rpcService.asPort().isOpen())
            rpcService.close();
        return true;
    }

public:
    /****************************************************************/
    Localizer() : closing(false) { }
};


/****************************************************************/
int main(int argc, char *argv[])
{
    Network yarp;
    ResourceFinder rf;
    rf.configure(argc,argv);

    if (!rf.check("file"))
    {
        if (!yarp.checkNetwork())
        {
            yError()<<"Unable to find Yarp server!";
            return EXIT_FAILURE;
        }
    }

    Localizer localizer;
    return localizer.runModule(rf);
}
