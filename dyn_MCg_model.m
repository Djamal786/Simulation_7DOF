function [M,C,g] = dyn_MCg_model(E_kin,E_pot,q,dq)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Mass matrix
M = sym(zeros(length(q)));
for n = 1:length(q)
    qn_coeffs = coeffs(E_kin,dq(n),'All');
    if length(qn_coeffs) == 3
        M(n,n) = 2*qn_coeffs(1);
        for m = n+1:length(q)
            M_nm = coeffs(qn_coeffs(2),dq(m),'All');
            if length(M_nm) == 2
                M(n,m) = M_nm(1);
                M(m,n) = M(n,m);
            end
        end
    end
end
% M = simplify(M);

% Coriolis matrix
C = sym(zeros(length(q)));
for i = 1:length(q)
    for j = 1:length(q)
        for k = 1:length(q)
            C(i,j) = C(i,j) + ( diff(M(i,j),q(k)) + diff(M(i,k),q(j)) - diff(M(k,j),q(i)) ) * dq(k);
        end
    end
end
C = 1/2 * C;%simplify(1/2 * C);

% Gravity
g = sym(zeros(size(q)));
for n = 1:length(q)
    g(n) = diff(E_pot, q(n));
end
% g = simplify(g);

end