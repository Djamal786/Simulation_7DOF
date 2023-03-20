function [omega] = cmp_omega(R,q,dq)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

R_dot = sym(zeros(size(R)));

for n = 1:length(q)
    R_dot = R_dot + diff(R,q(n)) * dq(n);
end

omega_tilde = R.' * R_dot;

omega = sym(zeros(3,1));
omega(1) = omega_tilde(3,2);
omega(2) = omega_tilde(1,3);
omega(3) = omega_tilde(2,1);

end