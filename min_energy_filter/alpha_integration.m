function A=alpha_integration(omega)

A=eye(3)-((1/2)*skew(omega))+(skew(omega)^2)/12;
end






