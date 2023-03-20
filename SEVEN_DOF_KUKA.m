clear
clc
close all

%%
coordinate_trafo_matrices;

%%
n = 7;
q = sym('q',[n 1]);
d = sym('d',[2 1]);
l = sym('l',size(q));
m = sym('m',size(q));
smplfy = 1; 

%% Forward kinematics of n-link
T_i = sym(zeros(4,4,length(q)));
for i = 1:length(q)
    T_i(:,:,i) = eye(4);
end
    T_i(:,:,1) = T_i(:,:,1) * rot_z(q(1));
    
    T_i(:,:,2) = T_i(:,:,2) * rot_z(q(1));
    T_i(:,:,2) = T_i(:,:,2) * rot_z(q(2)) * rot_x(sym(pi/2));  
    
    T_i(:,:,3) = T_i(:,:,3) * rot_z(q(1));
    T_i(:,:,3) = T_i(:,:,3) * rot_z(q(2))* rot_x(sym(pi/2));
    T_i(:,:,3) = T_i(:,:,3) * rot_z(q(3)) * rot_x(sym(-pi/2)); 

    T_i(:,:,4) = T_i(:,:,4) * rot_z(q(1));
    T_i(:,:,4) = T_i(:,:,4) * rot_z(q(2)) * rot_x(sym(pi/2));
    T_i(:,:,4) = T_i(:,:,4) * rot_z(q(3)) * rot_x(sym(-pi/2));
    T_i(:,:,4) = T_i(:,:,4) * rot_z(q(4)) * d_z(d(1)) * rot_x(sym(-pi/2));

    T_i(:,:,5) = T_i(:,:,5) * rot_z(q(1));
    T_i(:,:,5) = T_i(:,:,5) * rot_z(q(2)) * rot_x(sym(pi/2));
    T_i(:,:,5) = T_i(:,:,5) * rot_z(q(3)) * rot_x(sym(-pi/2));
    T_i(:,:,5) = T_i(:,:,5) * rot_z(q(4)) * d_z(d(1)) * rot_x(sym(-pi/2));
    T_i(:,:,5) = T_i(:,:,5) * rot_z(q(5)) * rot_x(sym(pi/2));

   
    T_i(:,:,6) = T_i(:,:,6) * rot_z(q(1));
    T_i(:,:,6) = T_i(:,:,6) * rot_z(q(2)) * rot_x(sym(pi/2));
    T_i(:,:,6) = T_i(:,:,6) * rot_z(q(3)) * rot_x(sym(-pi/2));
    T_i(:,:,6) = T_i(:,:,6) * rot_z(q(4)) * d_z(d(1)) * rot_x(sym(-pi/2));
    T_i(:,:,6) = T_i(:,:,6) * rot_z(q(5)) * rot_x(sym(pi/2));
    T_i(:,:,6) = T_i(:,:,6) * rot_z(q(6)) * d_z(d(2)) * rot_x(sym(pi/2));

    
    T_i(:,:,7) = T_i(:,:,7) * rot_z(q(1));
    T_i(:,:,7) = T_i(:,:,7) * rot_z(q(2)) * rot_x(sym(pi/2));
    T_i(:,:,7) = T_i(:,:,7) * rot_z(q(3)) * rot_x(sym(-pi/2));
    T_i(:,:,7) = T_i(:,:,7) * rot_z(q(4)) * d_z(d(1)) * rot_x(sym(-pi/2));
    T_i(:,:,7) = T_i(:,:,7) * rot_z(q(5)) * rot_x(sym(pi/2));
    T_i(:,:,7) = T_i(:,:,7) * rot_z(q(6)) * d_z(d(2)) * rot_x(sym(pi/2));
    T_i(:,:,7) = T_i(:,:,7) * rot_z(q(7)) * rot_x(sym(-pi/2));

    %Ist der Code zu aufwendig gestaltet oder gibt es irgend einen Fehler? 
    %Bis T_i(:,:,4) schafft mein Laptop es zu kompilieren, aber ab da
    %scheint er zu schwach zu sein. 

%% Velocities
dq = sym('dq',size(q));

% v = R dp/dt = R dp/dq dq/dt 
v_i = sym(zeros(3,length(q)));
for i = 1:length(q)
    v_i(:,i) = T_i(1:3,1:3,i) * jacobian(T_i(1:3,4,i),q) * dq;
end

% tilde(omega) = R^T dR/dt
omega_i = sym(zeros(3,length(q)));
for i = 1:length(q)
    omega_i(:,i) = cmp_omega(T_i(1:3,1:3,i),q,dq);
end

%% Kinetic energy

% link MOIs (alignment with principle axis coordinate system + rotational symmetries assumed)
I_i = sym(zeros(3,3,length(q)));
for i = 1:length(q)
    I_i(:,:,i) = diag(sym(['i_',num2str(i)],[3,1]));
end

% Kinetic energy (no eccentricities/displacements)
E_kin = sym(0);
for i = 1:length(q)
    E_kin = E_kin + (omega_i(:,i).')*I_i(:,:,i)*omega_i(:,i) + m(i)*(v_i(:,i).')*v_i(:,i);
end
E_kin = 1/2*E_kin;%simplify(1/2*E_kin);

%% Potential energy
% Gravitational potential (force in positive z-direction)
gp = sym('g');
U_g = 0;
for i = 1:length(q)
    l_sum = 0;
    for j = 1:i-1
        l_sum = l_sum + l(j);
    end
    l_sum = l_sum + l(i)/2;
    U_g = U_g + simplify(m(i) * gp * (l_sum - [0,0,1] * T_i(1:3,4,i)));
end
% m * g * (h_1 + h_2)
% h_1 = -x_1_fixed, x_1_fixed = x_B - r_1 --> h_1 = -x_B + r_1

%% Dynamics

[M,C,g] = dyn_MCg_model(E_kin,U_g,q,dq);

if smplfy
    M = simplify(M);
    C = simplify(C);
    g = simplify(g);
end
