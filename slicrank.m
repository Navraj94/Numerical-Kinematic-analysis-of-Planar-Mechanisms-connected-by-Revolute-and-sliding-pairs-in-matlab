function y=slicrank(param,x)
% a=[param(3)*cos(x(1)) -x(2); param(3)*sin(x(1)) 0];
% b=[-param(1)*cos(param(2));-param(1)*sin(param(2))-param(4)];
% y=a\b;
y=[param(1)*cos(param(2))+param(3)*cos(x(1))-x(2);param(1)*sin(param(2))+param(3)*sin(x(1))-param(4)];