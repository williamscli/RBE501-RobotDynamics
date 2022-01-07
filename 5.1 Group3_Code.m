% 3. Results
% 3.1 DH Parameters 
clear;clc;
syms theta1 theta2 theta3 theta4 theta5 theta6
% DH PARAMETERS
DHMatrix = [0,        90,  .089159, theta1;
            -.425,    0,   0,       theta2;
            -.39225,  0,   0,       theta3;
            0,        90,  .10915,  theta4;
            0,       -90,  .09465,  theta5;
            0,        0,   .0823,   theta6];
T06_DH=getTM(DHMatrix);
% 3.2 Forward Kinematics Using Product of Exponentials
theta_list=[theta1;theta2;theta3;theta4;theta5;theta6];
joint_transformations=cell(length(theta_list),1);

l1=0.089159;
l2=0.425;
l3=0.39225;
l4=0.10915;
l5=0.0946;
l6=0.0823;

I=eye(3);

v1=cross([1;0;0],[0;0;l1]);
v2=cross([1;0;0],[0;0;l1+l2]);
v3=cross([1;0;0],[0;0;l1+l2+l3]);
v4=cross([0;0;-1],[-l4;0;l1+l2+l3]);
v5=cross([1;0;0],[-l4;0;l1+l2+l3+l5]);
v6=cross([1;0;0],[-l4-l6;0;l1+l2+l3+l5]);

screw_axes = [[-1;0;0;v1],...
             [-1;0;0;v2],...
             [-1;0;0;v3],...
             [0;0;1;v4],...
             [-1;0;0;v5],...
             [-1;0;0;v6]];


for joint = 1 : length(theta_list)
    w = skew(screw_axes(1:3,joint));    
    v = screw_axes(4:6,joint);    
    theta = theta_list(joint);
    rotation = I + sin(theta)* w + (1-cos(theta))*w^2;    
    translation = (I * theta+(1-cos(theta))*w+(theta-sin(theta))*w^2)*v;    
    joint_transformations{joint} = [rotation translation; 0 0 0 1];
end

M06=[0 0 -1 (-l2-l3);
     -1 0 0 0;
     0 -1 0 (l1+l4+l5+l6);
     0 0 0 1];

T06_POE=simplify(joint_transformations{1}*joint_transformations{2}*joint_transformations{3}*joint_transformations{4}*joint_transformations{5}*joint_transformations{6}*M06);
% T=T06_POE-T06_DH
% 3.3 Jacobian

% JACOBIAN
JFinal = getJacobian(DHMatrix);

% Decoupling Jacobian
J11 = simplify(JFinal(1:3,1:3));
J22 = simplify(JFinal(4:6,4:6));

% Determinant calculations 
detJ11= simplify(det(J11));
detJ22= simplify(det(J22));
detJ = simplify(det(JFinal));

% 3.4 Dynamical Model
syms theta1(t) theta2(t) theta3(t) theta4(t) theta5(t) theta6(t) m1 m2 m3 m4 m5 m6 g t

%          Parameterize the Jacobian with time
JFinal     = subs(JFinal,[theta1 theta2 theta3 theta4 theta5 theta6],[theta1(t) theta2(t) theta3(t) theta4(t) theta5(t) theta6(t)]);

%          Define joint velocities
theta1_dot = diff(theta1(t),t);
theta2_dot = diff(theta2(t),t);
theta3_dot = diff(theta3(t),t);
theta4_dot = diff(theta4(t),t);
theta5_dot = diff(theta5(t),t);
theta6_dot = diff(theta6(t),t);

%          Homogenous transformation matrices
A1         = GetDh(DHMatrix(1,1), DHMatrix(1,2), DHMatrix(1,3), DHMatrix(1,4));
A2         = GetDh(DHMatrix(2,1), DHMatrix(2,2), DHMatrix(2,3), DHMatrix(2,4));
A3         = GetDh(DHMatrix(3,1), DHMatrix(3,2), DHMatrix(3,3), DHMatrix(3,4));
A4         = GetDh(DHMatrix(4,1), DHMatrix(4,2), DHMatrix(4,3), DHMatrix(4,4));
A5         = GetDh(DHMatrix(5,1), DHMatrix(5,2), DHMatrix(5,3), DHMatrix(5,4));
A6         = GetDh(DHMatrix(6,1), DHMatrix(6,2), DHMatrix(6,3), DHMatrix(6,4));

%          Linear velocities                        
v1         = JFinal(1:6,1)*theta1_dot;                                                                   v1sq       = v1.'*v1;
v2         = JFinal(1:6,1:2)*[theta1_dot,theta2_dot].';                                                  v2sq       = v2.'*v2;
v3         = JFinal(1:6,1:3)*[theta1_dot,theta2_dot,theta3_dot].';                                       v3sq       = v3.'*v3;
v4         = JFinal(1:6,1:4)*[theta1_dot,theta2_dot,theta3_dot,theta4_dot].';                            v4sq       = v4.'*v4;
v5         = JFinal(1:6,1:5)*[theta1_dot,theta2_dot,theta3_dot,theta4_dot,theta5_dot].';                 v5sq       = v5.'*v5;
v6         = JFinal(1:6,1:6)*[theta1_dot,theta2_dot,theta3_dot,theta4_dot,theta5_dot,theta6_dot].';      v6sq       = v6.'*v6;

%          Energies
k1         = 0.5.*m1.*v1sq;                                                                              p1         = m1.*g.*A1(3,4);
k2         = 0.5.*m2.*v2sq;                                                                              p2         = m2.*g.*A1(3,4)*A2(3,4);
k3         = 0.5.*m3.*v3sq;                                                                              p3         = m3.*g.*A1(3,4)*A2(3,4)*A3(3,4);
k4         = 0.5.*m4.*v4sq;                                                                              p4         = m4.*g.*A1(3,4)*A2(3,4)*A3(3,4)*A4(3,4);
k5         = 0.5.*m5.*v5sq;                                                                              p5         = m5.*g.*A1(3,4)*A2(3,4)*A3(3,4)*A4(3,4)*A5(3,4);
k6         = 0.5.*m6.*v6sq;                                                                              p6         = m6.*g.*A1(3,4)*A2(3,4)*A3(3,4)*A4(3,4)*A5(3,4)*A6(3,4);
k          = k1 + k2 + k3 + k4 + k5 + k6;                                                                p          = p1 + p2 + p3 + p4 + p5 + p6;

%          Lagrangian   
L          = k - p;

%          Torques
tau        = [diff(diff(L,theta1_dot),t) - diff(L,theta1(t));
              diff(diff(L,theta2_dot),t) - diff(L,theta2(t));
              diff(diff(L,theta3_dot),t) - diff(L,theta3(t));
              diff(diff(L,theta4_dot),t) - diff(L,theta4(t));
              diff(diff(L,theta5_dot),t) - diff(L,theta5(t));
              diff(diff(L,theta6_dot),t) - diff(L,theta6(t))];



% 3.5 Singularities
% SINGULARITIES
% Sweep for Elbow singularity
theta1 = 0; theta2 = -17;
theta3 = -135:5:70;
theta4 = 0; theta5 = 0.1; theta6 = 0;

X = theta3;
Y = eval(detJ);
SingPos=GetSing(X,Y);

plot(X,Y)
hold on
title('Elbow Singularity'); 
xlabel('\theta_3 Angle (degrees)'); 
ylabel('J_{11} Determinant');
plot(SingPos,0,'ro')
yline(0, 'r--')
text(SingPos(1),zeros(size(SingPos(1))),['   (', num2str(SingPos(1)), ', ', num2str(zeros(size(SingPos(1)))), ')'])

hold off

% Sweep for Shoulder singularity
theta1 = 5; theta2 = 30;
theta3 = -135:1:70;
theta4 = 0.1; theta5 = 0.1; theta6 = 5;

X = theta3;
Y = eval(detJ);
SingPos=GetSing(X,Y);

plot(X,Y)
title('Shoulder Singularity'); xlabel('\theta_3 Angle (degrees)'); ylabel('J Determinant');
hold on
xticks(-150:20:100)

theta2 = 60;
X2 = theta3;
Y2 = eval(detJ);
SingPos2=GetSing(X2,Y2);

plot(X2,Y2)
hold on

theta2 = -45;
X3 = theta3;
Y3 = eval(detJ);
SingPos3=GetSing(X3,Y3);

plot(X3,Y3)
plot(SingPos,0,'ro')
text(SingPos(1),zeros(size(SingPos(1))),['   (', num2str(SingPos(1)), ', ', num2str(zeros(size(SingPos(1)))), ')'])
plot(SingPos3(1),0,'ro')
text(SingPos3(1),zeros(size(SingPos3(1))),['   (', num2str(SingPos3(1)), ', ', num2str(zeros(size(SingPos3(1)))), ')'])
plot(SingPos2(2),0,'ro')
text(SingPos2(2),zeros(size(SingPos2(2))),['   (', num2str(SingPos2(2)), ', ', num2str(zeros(size(SingPos2(2)))), ')'])
yline(0, 'r--')
legend('\theta_2 = 30^{\circ}', '\theta_2 = 60^{\circ}', '\theta_2 = -45^{\circ}')
hold off

% Sweep for wrist singularity
theta1 = 30; theta2 = 30; theta3 = 30;
theta4 = 30; theta6 = 30;
theta5 = 0:5:360;

X=theta5;
Y=eval(detJ22);
SingPos=GetSing(X,Y); % Get Singularity positions

plot(X,Y)
hold on
title('Wrist Singularity'); 
xlabel('\theta_5 Angle (degrees)'); 
ylabel('J_{22} Determinant');
plot(SingPos,0,'ro')
text(SingPos(1),zeros(size(SingPos(1))),['   (', num2str(SingPos(1)), ', ', num2str(zeros(size(SingPos(1)))), ')'])
text(SingPos(2),zeros(size(SingPos(2))),['   (', num2str(SingPos(2)), ', ', num2str(zeros(size(SingPos(2)))), ')'])
text(SingPos(3),zeros(size(SingPos(3))),['   (', num2str(SingPos(3)), ', ', num2str(zeros(size(SingPos(3)))), ')'])
yline(0, 'r--')
hold off
% 3.6 Motion Planning
%Trajectory/Path Planning
%Importing URDF file for the selected Universal Robot with designed EOA
%Imported robot would be named kitchenRobot
kitchenRobot = importrobot('UR5_fliper.urdf');

%InteractiveRigidBodyTree creates a fig that shows the imported robot and
%enables us to interactively modify the robot configuration using an
%interactive marker
iviz_robot = interactiveRigidBodyTree(kitchenRobot);
iviz_robot.MarkerBodyName = "L6"; %would be used to determine joint configs
iviz_robot.ShowMarker = false;   

%%%%%%%%%%%%%% Environment setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = gca;                                                     %%%%%%%%
xlim([-2 2])                                                  %%%%%%%%
ylim([-2 2])                                                  %%%%%%%%
zlim([-1 1])                                                  %%%%%%%%
                                                              %%%%%%%%
% Table platform                                              %%%%%%%%
platform = collisionBox(2.5,3,1);                             %%%%%%%%
platform.Pose = trvec2tform([0 0 -0.5]);                      %%%%%%%%
show(platform,'Parent', ax);                                  %%%%%%%%
                                                              %%%%%%%%
%Left surafce                                                 %%%%%%%%
leftSurface = collisionBox(1,0.5,0.05);                       %%%%%%%%
leftSurface.Pose = trvec2tform([-0.05 -0.5 0.025]);           %%%%%%%%
[~, patchObj] = show(leftSurface,'Parent',ax);                %%%%%%%%
patchObj.FaceColor = [1 1 1];                                 %%%%%%%%
                                                              %%%%%%%%
%Right surface where raw product is picked.                   %%%%%%%%
rightSurface = collisionBox(1,0.5,0.05);                      %%%%%%%%
rightSurface.Pose = trvec2tform([-0.05 0.75 0.025]);          %%%%%%%%
[~, patchObj] = show(rightSurface,'Parent',ax);               %%%%%%%%
patchObj.FaceColor = [1 1 1];                                 %%%%%%%%
                                                              %%%%%%%%
%Grilling Surface                                             %%%%%%%%
grillSurface = collisionBox(0.5,1,0.05);                      %%%%%%%%
grillSurface.Pose = trvec2tform([0.75 0 0.025]);              %%%%%%%%
[~, patchObj] = show(grillSurface,'Parent',ax);               %%%%%%%%
patchObj.FaceColor = [1 0 0];                                 %%%%%%%%
                                                              %%%%%%%%
                                                              %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Trajectory Generation using desired joint configurations and performing
%inverse Kinematics


 jointConfigurations = [-2.47001158713855,-2.47166292992973,-0.936323863734110,...
     -0.936323863734110,-0.936323863734110,0.583352843504152,0.583352843504152,...
     2.32433046879710; 1.43537931847073,1.75626278874233,2.26006376954913,...
     2.26006376954913,2.26006376954913,1.34740311212141,1.34740311212141,...
     1.34740311212141; -0.127576544429983,-0.445549605138962,-1.51404359064025,...
     -1.51404359064025,-1.51404359064025,-0.165141294576372,-0.165141294576372,...
     0.175563603955479; -1.13656830415248,-1.21439452119498,-0.670632287531426,...
     -0.670632287531426,-0.670632287531426,-1.09849201643661,-1.09849201643661,...
     -3.19589355327972; 0.00281931715369363,0.000281821285244265,0.000281821285244265,...
     0.000281821285244265,0.000281821285244265,0.000281821285244265,0.000281821285244265,...
     0.000281821285244265; 0,0,0,3.11839689699530,6.25272436852864,...
     6.29083090834297,7.46876093645159,6.24190791895185];

% A trapezoidal velocity profile is generated using trapveltraj. 
% means the robot stops smoothly at each waypoint, 
% but achieves a set max speed while in motion.

waypointsSamples = 100*size(jointConfigurations, 2) + 1;
[q,~,~, tvec] = trapveltraj(jointConfigurations,waypointsSamples,'EndTime',10);

showFigure(iviz_robot)
%Replay the generated trajectory by iterating the generated q matrix, 
%which represents a series of joint configurations that move between each waypoint
rateCtrlObj = rateControl(waypointsSamples/(max(tvec) + tvec(2)));
for i = 1:waypointsSamples
    iviz_robot.Configuration = q(:,i);
    waitfor(rateCtrlObj);
end



% 6. Appendix/Functions
% This function computes the homgeneous transformation matrix.
function A = GetDh(a,alpha,d,theta)
%              a               aplha                        d                    theta
         A=[cosd(theta)  -sind(theta)*cosd(alpha)  sind(theta)*sind(alpha)   a*cosd(theta);
            sind(theta)  cosd(theta)*cosd(alpha)   -cosd(theta)*sind(alpha)  a*sind(theta);
            0            sind(alpha)               cosd(alpha)               d;
            0            0                         0                         1];
end

% This function takes in a 6x4 DH parameters matrix for a robot with six revolute joints and outputs the resulting transformation matrix of the end effector.
function T06 = getTM(DHMatrix)
% Function for calculating transformation matrix for EOF with 6 revolute joints using a DH table
% transformation matrices calculations
syms theta1 theta2 theta3 theta4 theta5 theta6
A1 = GetDh(DHMatrix(1,1), DHMatrix(1,2), DHMatrix(1,3), DHMatrix(1,4));
A2 = GetDh(DHMatrix(2,1), DHMatrix(2,2), DHMatrix(2,3), DHMatrix(2,4));
A3 = GetDh(DHMatrix(3,1), DHMatrix(3,2), DHMatrix(3,3), DHMatrix(3,4));
A4 = GetDh(DHMatrix(4,1), DHMatrix(4,2), DHMatrix(4,3), DHMatrix(4,4));
A5 = GetDh(DHMatrix(5,1), DHMatrix(5,2), DHMatrix(5,3), DHMatrix(5,4));
A6 = GetDh(DHMatrix(6,1), DHMatrix(6,2), DHMatrix(6,3), DHMatrix(6,4));
T06 = simplify(A1*A2*A3*A4*A5*A6);
end
% This function takes in a 6x4 DH parameters matrix for a robot with six revolute joints and outputs the resulting Jacobian of the robot.
function J = getJacobian(DHMatrix)
% Function for calculating Jacobian for a robot with 6 revolute joints using a DH table
% transformation matrices calculations
syms theta1 theta2 theta3 theta4 theta5 theta6
A1 = GetDh(DHMatrix(1,1), DHMatrix(1,2), DHMatrix(1,3), DHMatrix(1,4));
A2 = GetDh(DHMatrix(2,1), DHMatrix(2,2), DHMatrix(2,3), DHMatrix(2,4));
A3 = GetDh(DHMatrix(3,1), DHMatrix(3,2), DHMatrix(3,3), DHMatrix(3,4));
A4 = GetDh(DHMatrix(4,1), DHMatrix(4,2), DHMatrix(4,3), DHMatrix(4,4));
A5 = GetDh(DHMatrix(5,1), DHMatrix(5,2), DHMatrix(5,3), DHMatrix(5,4));
A6 = GetDh(DHMatrix(6,1), DHMatrix(6,2), DHMatrix(6,3), DHMatrix(6,4));
T02 = simplify(A1*A2);
T03 = simplify(A1*A2*A3);
T04 = simplify(A1*A2*A3*A4);
T05 = simplify(A1*A2*A3*A4*A5);
T06 = simplify(A1*A2*A3*A4*A5*A6);

% Origin calculations
o0 = [0; 0; 0];
o1 = A1(1:3, 4);
o2 = T02(1:3, 4);
o3 = T03(1:3, 4);
o4 = T04(1:3, 4);
o5 = T05(1:3, 4);
o6 = T06(1:3, 4);

% z calculations
z0 = [0; 0; 1];
z1 = A1(1:3, 3);
z2 = T02(1:3, 3);
z3 = T03(1:3, 3);
z4 = T04(1:3, 3);
z5 = T05(1:3, 3);
z6 = T06(1:3, 3);

% Jacobian calculations
J1 = [cross(z0, (o6-o0)); z0];
J2 = [cross(z1, (o6-o1)); z1];
J3 = [cross(z2, (o6-o2)); z2];
J4 = [cross(z3, (o6-o3)); z3];
J5 = [cross(z4, (o6-o4)); z4];
J6 = [cross(z5, (o6-o5)); z5];
J = [J1 J2 J3 J4 J5 J6];
end

% This function takes theta array (X) and detJ as Y and returns the theta values where singularity occurs. 
function SingularityPos=GetSing(X,Y)
% Find thetas where Jacobian is 0
SingularityPos=[];
YRound=round(Y,6);
for i=1:length(YRound)
    if YRound(:,i)==0.00000
        SingularityPos(:,end+1)=X(:,i);
    end
end
disp("Signularities occur at theta(s): " + num2str(SingularityPos))
end