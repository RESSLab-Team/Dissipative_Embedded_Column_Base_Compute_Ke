clear all; clc; close all;

%% Read me
% This code calculates the elastic stiffness of dissipative embedded column
% base connections for given input parameters. The users are referred to the references [1] (chapter 5 and appendix) and [2] for derivation of the equations. 
%
% Written by Hiroyuki Inamasu, (RESSLab, Ecole Polytechnique federale de Lausanne, EPFL)  
%
% Last modify: 12 April 2021 
%
%
% References:
% 
% [1] H Inamasu (2021), Development of dissipative steel column-base
% connections for steel moment-resisting frames under seismic loading, EPFL
% PhD thesis.
%
% [2] H Inamasu, A de Castro e Sousa, DG Lignos (2021), Development and experimental
% validation of dissipative embedded column base connections for enhanced seismic performance of steel moment resisting frames
% ASCE Journal of structural engineering (under review).


%% Inputs
% Define column cross section
% Dimensions of unreduced section
Col_d = 400.0;% mm. Depth 
Col_bf = 180.0;% mm. Flange width
Col_tw = 8.6;% mm. Web thickness
Col_tf = 13.5;% mm. Flange thickness
Col_r_corner = 21.0;% mm. Raidus of k-area
% Dimensions of reduced section
Col_d_reduced = 400.0;% mm. Depth
Col_bf_reduced = 88.0;% mm. Flange width
Col_tw_reduced = 8.6;% mm. Web thickness
Col_tf_reduced = 13.5;% mm. Flange thickness
Col_r_corner_reduced = 21.0;% mm. Raidus of k-area

% Define steel material properties
Col_E = 200000; % MPa. Young's modulus of the column steel material
Col_G = Col_E/2/(1+0.3); % MPa. Shear elastic modulus of the column steel material

% Define geometric parameters 
L_c = 1525; % mm. Inflection point height of the column measured from the foundation top surface to the column inflection point 
d_embed = 850; % mm. Embedment depth of the column
e_dt = 100; % mm. Top unreduced length of the embedded column portion 
r_cut = 46; % mm. Radius of the flange cut with quadrant shape
H = L_c+d_embed; % mm. Inflection point height + Embedment depth
L = d_embed-e_dt/2; % mm. Distance from the column inflection point to the level of bearing reaction forces

% Define debonding material properties
t_rubber = 1; % mm. Thickness of the debonding material
Num_BearingSurface = 2; % Number of bearing surfaces
Bearing_A = Col_bf*e_dt; % mm2. Effective bearing area
E_rubber = 1; % MPa 
k_rubber = E_rubber*Bearing_A/t_rubber*Num_BearingSurface; % N/mm


% Compute column section properties
% Unreduced section
Col_A  = 2*Col_bf*Col_tf+Col_tw*(Col_d-2*Col_tf)+0.8584*Col_r_corner^2; % Cross-sectional area
Col_I  = (Col_bf*Col_d^3-(Col_bf-Col_tw)*(Col_d-2*Col_tf)^3)/12.0+0.8584*Col_r_corner^2*(0.5*Col_d-Col_tf-0.4467*Col_r_corner/2)^2;% Moment of inertia about strong axis
Col_A_v = Col_A-2*Col_bf*Col_tf+(Col_tw+2*Col_r_corner)*Col_tf;% shear area
% Reduced section
Col_I_reduced  = (Col_bf_reduced*Col_d_reduced^3-(Col_bf_reduced-Col_tw_reduced)*(Col_d_reduced-2*Col_tf_reduced)^3)/12.0+0.8584*Col_r_corner_reduced^2*(0.5*Col_d_reduced-Col_tf_reduced-0.4467*Col_r_corner_reduced/2)^2;% Moment of inertia about strong axis



%% Compute elastic stiffness of dissipative embedded column base connections
% Assume unit base moment applied to the cantilever column
M_kNm = 1; % kNm
M = M_kNm*1000*1000; % Nmm
P = M/L_c; % N

Error_min = 10^1000;
Error=zeros(1,1000);

% Set distances for integration
d_1 = e_dt+(r_cut/2);
d_2 = e_dt/2+(r_cut/2);

Max_x_i=10000; Squared_Error = zeros(1,Max_x_i);
for x_i = 1:Max_x_i

    % Assume normal displacement of debonding material
    x_rubber = x_i/10000; % mm.   
    

    % Compute lateral deflection and rotation at the level of bearing
    % reaction forces under an applied unit base moment through integration
    
    % Boundary conditions at the level of the base plate
    % Boundary condition at y_depth(1) = 0
    y_depth(1) = 0;
    dv_dy(1) = 0;
    v(1)     = 0;

    % C_1 and C_2 for range y_depth(1) = 0 < y < y_depth(2) = d_1
    I_now = Col_I;
    C_1(1) = dv_dy(1)-(1/Col_E/I_now*(1/2*(k_rubber*x_rubber-P)*y_depth(1)^2+(P*H-k_rubber*x_rubber*L)*y_depth(1))); 
    C_2(1) = v(1)-(1/Col_E/I_now*(1/6*(k_rubber*x_rubber-P)*y_depth(1)^3+1/2*(P*H-k_rubber*x_rubber*L)*y_depth(1)^2))-C_1(1)*y_depth(1); 

    % Boundary condition at y_depth(2) = d_1
    y_depth(2) = d_1;
    dv_dy(2) = 1/Col_E/I_now*(1/2*(k_rubber*x_rubber-P)*y_depth(2)^2+(P*H-k_rubber*x_rubber*L)*y_depth(2))+C_1(1);
    v(2)     = 1/Col_E/I_now*(1/6*(k_rubber*x_rubber-P)*y_depth(2)^3+1/2*(P*H-k_rubber*x_rubber*L)*y_depth(2)^2)+C_1(1)*y_depth(2)+C_2(1);

    % C_1 and C_2 for range y_depth(2) = d_1 < y < y_depth(3) = L-d_2/2
    I_now = Col_I_reduced;
    C_1(2) = dv_dy(2)-(1/Col_E/I_now*(1/2*(k_rubber*x_rubber-P)*y_depth(2)^2+(P*H-k_rubber*x_rubber*L)*y_depth(2))); 
    C_2(2) = v(2)-(1/Col_E/I_now*(1/6*(k_rubber*x_rubber-P)*y_depth(2)^3+1/2*(P*H-k_rubber*x_rubber*L)*y_depth(2)^2))-C_1(2)*y_depth(2); 

    % Boundary condition at y_depth(3) = L-d_2/2
	y_depth(3) = L-d_2/2;
    dv_dy(3) = 1/Col_E/I_now*(1/2*(k_rubber*x_rubber-P)*y_depth(3)^2+(P*H-k_rubber*x_rubber*L)*y_depth(3))+C_1(2);
    v(3)     = 1/Col_E/I_now*(1/6*(k_rubber*x_rubber-P)*y_depth(3)^3+1/2*(P*H-k_rubber*x_rubber*L)*y_depth(3)^2)+C_1(2)*y_depth(3)+C_2(2);

    % C_1 and C_2 for range y_depth(3) = L-d_2/2 < y < y_depth(4) = L
    I_now = Col_I;
    C_1(3) = dv_dy(3)-(1/Col_E/I_now*(1/2*(k_rubber*x_rubber-P)*y_depth(3)^2+(P*H-k_rubber*x_rubber*L)*y_depth(3))); 
    C_2(3) = v(3)-(1/Col_E/I_now*(1/6*(k_rubber*x_rubber-P)*y_depth(3)^3+1/2*(P*H-k_rubber*x_rubber*L)*y_depth(3)^2))-C_1(3)*y_depth(3); 

    % Boundary condition at y_depth(4) = L
	y_depth(4) = L;
    I_now = Col_I;
    dv_dy(4) = 1/Col_E/I_now*(1/2*(k_rubber*x_rubber-P)*y_depth(4)^2+(P*H-k_rubber*x_rubber*L)*y_depth(4))+C_1(3);
    v(4)     = 1/Col_E/I_now*(1/6*(k_rubber*x_rubber-P)*y_depth(4)^3+1/2*(P*H-k_rubber*x_rubber*L)*y_depth(4)^2)+C_1(3)*y_depth(4)+C_2(3);


    % Lateral deflection due to flexure at the level of bearing reaction forces 
    Delta_f_gl = v(4); % mm.
    % Rotation due to flexure at the level of bearing reaction forces 
    Theta_gl = dv_dy(4); % rad.
    % Lateral deflection due to shear at the level of bearing reaction forces 
    Delta_s_gl = (P-k_rubber*x_rubber)/(Col_G*Col_A_v/L);% mm
    
    % Calculate squared error between assumed normal displacement of debonding material and normal displacement of debonding material under an applied unit base moment   
    Squared_Error(x_i) = (x_rubber-(Delta_f_gl+Delta_s_gl))^2;
    
    % Obtain lateral deflection and rotation that provide minimum error 
    if Squared_Error(x_i) < Error_min 
        Error_min = Squared_Error(x_i);
        x_last = x_rubber;
        Theta_last = Theta_gl;
    end
    

end

% Total lateral deflection of the cantilever column at the column top under an
% applied base moment
Delta_at_P = P/(3*Col_E*Col_I/(H-L)^3)+P/(Col_G*Col_A_v/(H-L))+x_last+Theta_last*(H-L); % mm

% Column top lateral deflection due to dissipative embedded column base flexibility under an
% applied base moment (Total lateral deflection - column lateral deformation)
Delta_ar_P_due_to_DECB = Delta_at_P-P/(3*Col_E*Col_I/(L_c)^3)-P/(Col_G*Col_A_v/(L_c));

% Total stiffness in lateral force - column top lateral displacement relation 
K_V_dh = P/Delta_at_P; % N/mm. (V-delta_horizontal)

% Total stiffness in base moment - column drift ratio relation 
K_M_theta = K_V_dh*L_c^2; % Nmm/rad. (M_base-theta)

% Total stiffness in base moment - column drift ratio relation 
K_M_theta_kNm_rad = K_M_theta/1000/1000; % kNm/rad. (M_base-theta)

% Stiffness of dissipative embedded base (DECB) in base moment - DECB rotation relation 
K_M_theta_DECBspring_kNm_perrad = (P/Delta_ar_P_due_to_DECB*L_c^2)/1000/1000; % kNm/rad. 


%% Print information to command window
fprintf('Results: %g\n','')
fprintf(' %g\n','')
fprintf('Total stiffness in base moment - column drift ratio relation: %g (kNm/rad) \n',K_M_theta_kNm_rad)
fprintf(' %g\n','')
fprintf('Rotational stiffness of dissipative embedded base (DECB) (stiffness in base moment - DECB rotation relation): %g (kNm/rad) \n',K_M_theta_DECBspring_kNm_perrad)
fprintf(' %g\n','')

