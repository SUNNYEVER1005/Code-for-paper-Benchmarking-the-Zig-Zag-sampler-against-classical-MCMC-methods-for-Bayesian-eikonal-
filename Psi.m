% potential in Zig-Zag
function value=Psi(pos,timeobs,sigma_r,receiver,v0)

T_value=sqrt( ( pos(1)-receiver(:,1) ).^2+( pos(2)-receiver(:,2)).^2 )/v0;

value=0.5*sum(   (T_value-timeobs+pos(3)).^2 ) /sigma_r^2;

end