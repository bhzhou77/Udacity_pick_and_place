## Project: Kinematics Pick & Place
For more information, refer to [this repository](https://github.com/udacity/RoboND-Kinematics-Project).

[//]: # (Image References)

[image1]: figure_kr210.jpg
[image2]: homogeneous_transformation.jpg
[image3]: individual_homo_trans.jpg
[image4]: total_homo_trans_euler_angle.jpg
[image5]: wrist_center.jpg
[image6]: first_three_geometry.jpg
[image7]: first_three_derivation.jpg
[image8]: first_three_rotation.jpg
[image9]: last_three_derivation.jpg

### Kinematic Analysis
#### 1. Run the forward_kinematics demo and evaluate the kr210.urdf.xacro file to perform kinematic analysis of Kuka KR210 robot and derive its DH parameters.

The illustrative figure of the kr210 arm is shown below in Fig. 1, where the
assignment of link and joint numbers, the origins and orientations of the 
reference frames, and the parameters are consistent with the modified DH
convention (Craig's, 2005).

![Figure 1][image1]

With the definitions of DH parameters in Fig. 1 and informations shown in the
urdf file, I obtained the following DH parameter table.

Links | alpha(i-1) | a(i-1) | d(i) | theta(i)
--- | --- | --- | --- | ---
0->1 | 0 | 0 | 0.75 | q1
1->2 | - pi/2 | 0.35 | 0 | q2 - pi/2
2->3 | 0 | 1.25 | 0 | q3
3->4 |  - pi/2 | - 0.054 | 1.5 | q4
4->5 | pi/2 | 0 | 0 | q5
5->6 | - pi/2 | 0 | 0 | q6
6->EE | 0 | 0 | 0.303 | 0

#### 2. Using the DH parameter table you derived earlier, create individual transformation matrices about each joint. In addition, also generate a generalized homogeneous transform between base_link and gripper_link using only end-effector(gripper) pose.

In general, a homogeneous transformation matrix is shown in Fig. 2, which can
generate the individual transformation matrix for any joint (Fig. 3) if being 
fed with corresponding row of numbers in the DH table.

![Figure 2][image2]

![Figure 3][image3]

Under Tait-Bryan convention, the yaw, pitch, and roll are defined as the
intrinsic rotation Z1-Y2-X3, and they can be labelled as alpha, beta, and gamma.
This intrinsic rotation is done with respect to the reference frame of the end 
effector under the URDF definition, the initial orientation of which is the
same as the base link. Thus, the homogeneous transformation matrix from the 
base link to the gripper link using only the position and orientation of the
gripper link has the form shown in Fig. 4, where alpha, beta and gamma are the
yaw, pitch and roll values from ros, and px, py, and pz are the three position
components (in base link reference frame) from ros.

![Figure 4][image4]

#### 3. Decouple Inverse Kinematics problem into Inverse Position Kinematics and inverse Orientation Kinematics; doing so derive the equations to calculate all individual joint angles.

The labelled R_total in Fig. 4 represents the total rotation of the 
URDF-defined gripper reference frame with respect to the base link reference
frame. I will do two body-fixed intrinsic rotation on R_total: the first one is
180 degree on Z1, represented by R_Z(180), and the second is -90 degree on Y2,
represented by R_Y(-90). After this, the displacement of the end effector d_G is
on the Z_6 or Z_G direction. Thus, I can get the position of the wrist center in
the base link frame as shown in Fig. 5.

![Figure 5][image5]

Once I have the position of the wrist center, I can derive the first three joint
angles (Fig. 6 and Fig. 7).

![Figure 6][image6]

![Figure 7][image7]

Then, I extract the first three rotation matrices from Fig. 3, and evaluate 
their product using the obtained theta1, theta2 and theta 3 (Fig. 8).

![Figure 8][image8]

With this, I can get the total rotation from joint 3 to joint 6, and derive
the last three joint angles (Fig. 9).

![Figure 9][image9]



### Project Implementation

#### 1. Fill in the `IK_server.py` file with properly commented python code for calculating Inverse Kinematics based on previously performed Kinematic Analysis. Your code must guide the robot to successfully complete 8/10 pick and place cycles. Briefly discuss the code you implemented and your results. 

First, I define all the DH related parameters and matrices with respect to the
modified DH convention (Craig's, 2005) as shown below.

```{python}
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

# Total rotational matrix from base link to link 3
R0_3 = T0_1[0:3,0:3] * T1_2[0:3,0:3] * T2_3[0:3,0:3]
```

Then, I defined the corrected total rotation matrix of EE:

```{python}
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
```

After all these definitions, I define those helper functions that I am going to
use to solve the IK problem.

```{python}
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
```

At last, inside the function `handle_calculate_IK`, I calculate the six joint
angles given the position and orientation of EE.

```{python}
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
```

I tried this code several times, and the arm could always successfully complete 
the pick-and-place task. I defintely need a better computer for this, since with
my current one, the arm moves in a slow motion and not very smooth sometimes.

