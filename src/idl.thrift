/******************************************************************************
 *                                                                            *
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file idl.thirft
 * @authors: Giulia Vezzani <giulia.vezzani@iit.it>
 */

 /**
 * Vector
 *
 * IDL structure to set/show advanced parameters.
 */
 struct Vector
 {
 } (
    yarp.name = "yarp::sig::Vector"
    yarp.includefile="yarp/sig/Vector.h"
 )

 /**
 * Bottle
 *
 * IDL structure to set/show advanced parameters.
 */
 struct Bottle
 {
 } (
    yarp.name = "yarp::os::Bottle"
    yarp.includefile="yarp/os/Bottle.h"
 )

 /**
 * PointCloudXYZRGBA
 *
 * IDL structure to set/show advanced parameters.
 */
 struct PointCloudXYZRGBA
 {
 }
	(
    yarp.name = "yarp::sig::PointCloudXYZRGBA"
    yarp.includefile="yarp/sig/PointCloud.h"
 )

 service Localizer_IDL
 {
	 /**
     * Localize the superquadric pose (position and Euler angles)
     * @param points is a list of DataXYZ of the object
     * @param object_name is the object name
     * @return a 5D Vector containing the dimensions and shape
     */
	 Vector localize_superq(1: string &object_names, 2: Bottle &points, 3: i32 object_num);

     /**
     * Set the options to remove outliers from point cloud
     * @param values is a bottle containing the number of points
     * and the radius used by the filter
     * @return true
     */
	 bool set_remove_outliers(1: Bottle &values);

	 /**
     * Set the option to uniformly sample the point cloud on/off
     * @param input can be on or off
     * @return true
     */
     bool set_uniform_sample(2: string &input);

	 /**
     * Set the option to randomly sample the point cloud on/off
     * @param input can be on or off
     * @return true
     */
	 bool set_random_sample(2: string &input);

	 /**
     * Set the option analytic_gradient for the optimization problem on/off
     * @param input can be on or off
     * @return true
     */
	 bool set_analytic_gradient(2: string &input);

     /**
     * Reset the visualization
     * @return true
     */
     bool reset_visualization();

 }
