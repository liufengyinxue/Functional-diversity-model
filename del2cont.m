function d2x=del2cont(x,n) % receive a vector and return a vector of 2nd "spatial" derivative.
% Last chack: 14.6.09
if n<2
    d2x(1:n)=0;
    return
end
d2x(1)=(x(2)-x(1))*n^2;
d2x(2:n-1)=(x(1:n-2)-2*x(2:n-1)+x(3:n))*n^2;
d2x(n)=(x(n-1)-x(n))*n^2;

return;
