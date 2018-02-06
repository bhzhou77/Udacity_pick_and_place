#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *

##################################################################
############### Modified DH parameters and matrices ##############
##################################################################

# Define DH parameter symbols
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7') # twist angle
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7') # link length
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8') # link offset
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') # joint angle

# Modified DH parameters
DH_table = {alpha0:      0, a0:      0, d1:   0.75, q1:        q1,
            alpha1: -pi/2., a1:   0.35, d2:      0, q2: -pi/2.+q2,
            alpha2:      0, a2:   1.25, d3:      0, q3:        q3,
            alpha3: -pi/2., a3: -0.054, d4:    1.5, q4:        q4,
            alpha4:  pi/2., a4:      0, d5:      0, q5:        q5,
            alpha5: -pi/2., a5:      0, d6:      0, q6:        q6,
            alpha6:      0, a6:      0, d7:  0.303, q7:         0}

# Modified DH transformation matrix
def TF_matrix(alpha, a, d, q):
    TF = Matrix([[           cos(q),           -sin(q),           0,             a],
                 [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                 [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                 [                0,                 0,           0,             1]])
    return TF

# Individual transformation matrices
T0_1  = TF_matrix(alpha0, a0, d1, q1).subs(DH_table)
T1_2  = TF_matrix(alpha1, a1, d2, q2).subs(DH_table)
T2_3  = TF_matrix(alpha2, a2, d3, q3).subs(DH_table)
T3_4  = TF_matrix(alpha3, a3, d4, q4).subs(DH_table)
T4_5  = TF_matrix(alpha4, a4, d5, q5).subs(DH_table)
T5_6  = TF_matrix(alpha5, a5, d6, q6).subs(DH_table)
T6_EE = TF_matrix(alpha6, a6, d7, q7).subs(DH_table)

# Total transformation matrix from base link to end effector
T0_EE = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_EE

# Total rotational matrix from base link to link 3
R0_3 = T0_1[0:3,0:3] * T1_2[0:3,0:3] * T2_3[0:3,0:3]


###############################################################
############# Corrected EE rotation matrix ####################
###############################################################

r, p, y = symbols('r p y') 

ROT_x = Matrix([[1,      0,       0],
                [0, cos(r), -sin(r)],
                [0, sin(r), cos(r)]])

ROT_y = Matrix([[ cos(p), 0, sin(p)],
                [      0, 1,      0],
                [-sin(p), 0, cos(p)]])

ROT_z = Matrix([[cos(y), -sin(y), 0],
                [sin(y),  cos(y), 0],
                [     0,       0, 1]])

ROT_EE = ROT_z * ROT_y * ROT_x

# Account for the URDF file defintion and make the correction.
ROT_correction = ROT_z.subs(y, radians(180)) * ROT_y.subs(p, radians(-90))
ROT_EE = ROT_EE * ROT_correction


###############################################################
############# Helper functions ################################
###############################################################

def get_corrected_EE_rotation(roll, pitch, yaw):
    return ROT_EE.subs({r: roll, p: pitch, y: yaw})


def get_WC(px, py, pz, ROT_EE):
    EE = Matrix([[px],
                 [py],
                 [pz]])

    WC = EE - 0.303 * ROT_EE[:, 2]

    return WC


def get_first_three_joint_angles(WC):

    # Joint angles through IK
    theta1 = atan2(WC[1], WC[0])

    # SSS triangle for theta2 and theta3
    side_a = sqrt(0.054**2+1.5**2)
    side_b = sqrt((sqrt(WC[1]**2+WC[0]**2)-0.35)**2+(WC[2]-0.75)**2)
    side_c = 1.25

    angle_a = acos((side_b**2+side_c**2-side_a**2)/(2*side_b*side_c))
    angle_b = acos((side_a**2+side_c**2-side_b**2)/(2*side_a*side_c))
    angle_c = acos((side_a**2+side_b**2-side_c**2)/(2*side_a*side_b))

    theta2 = pi/2. - angle_a - atan2(WC[2]-0.75,sqrt(WC[0]**2+WC[1]**2)-0.35)
    theta3 = pi/2. - angle_b - atan2(0.054,1.5)

    return theta1, theta2, theta3


def evaluate_R0_3(theta1, theta2, theta3):
    return R0_3.evalf(subs={q1:theta1, q2:theta2, q3:theta3})


def get_last_three_joint_angles(theta1, theta2, theta3, R0_3, ROT_EE):

    R3_6 = R0_3.T * ROT_EE

    # Euler angles from rotation matrix
    theta4 = atan2(R3_6[2,2], -R3_6[0,2])
    theta5 = atan2(sqrt(R3_6[0,2]**2+R3_6[2,2]**2), R3_6[1,2])
    theta6 = atan2(-R3_6[1,1], R3_6[1,0])

    return theta4, theta5, theta6



def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            # Corrected EE rotation matrix
            ROT_EE = get_corrected_EE_rotation(roll, pitch, yaw)
            # Wrist center
            WC = get_WC(px, py, pz, ROT_EE)
            # First three joint angles
            theta1, theta2, theta3 = get_first_three_joint_angles(WC)
            # Evaluate R0_3
            R0_3 = evaluate_R0_3(theta1, theta2, theta3)
            # Last three joint angles
            theta4, theta5, theta6 = get_last_three_joint_angles(theta1, theta2, theta3, R0_3, ROT_EE)
            

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)
    
        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
