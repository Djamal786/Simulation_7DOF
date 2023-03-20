syms r_x(theta) r_y(theta) r_z(theta) rot_x(theta) rot_y(theta) rot_z(theta) tran(d) d_x(delta)

r_x(theta) = [1,0,0;0,cos(theta),-sin(theta);0,sin(theta),cos(theta)];
r_y(theta) = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
r_z(theta) = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];

rot_x(theta) = [r_x(theta), zeros(3,1); zeros(1,3), 1];
rot_y(theta) = [r_y(theta), zeros(3,1); zeros(1,3), 1];
rot_z(theta) = [r_z(theta), zeros(3,1); zeros(1,3), 1];



d = sym('d',[3 1]);
tran(d) = [eye(3), d;  zeros(1,3), 1];

d_x(delta) = tran(delta,0,0);
d_y(delta) = tran(0,delta,0);
d_z(delta) = tran(0,0,delta);


