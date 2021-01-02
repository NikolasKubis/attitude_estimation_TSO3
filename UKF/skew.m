function tildematrix=skew(m)

% The function skew performs an isomorphism between the 3x1 arrays and the
% 3x3 skew symmetric matrices .

m1=m(1);
m2=m(2); 
m3=m(3);
tildematrix=[0 -m3 m2;m3 0 -m1;-m2 m1 0];
end