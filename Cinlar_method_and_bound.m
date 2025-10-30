%calculation of the M and tau
function [M_value,tau]=Cinlar_method_and_bound(vel,t,pos,timeobs,sigma_r,gamma,receiver,v0)
dim=size(pos);
u=rand(dim);
s=-log(u);
a=zeros(dim(1),1);b=a;

%local maximum value
T_value=sqrt( ( pos(1)+vel(1)*t-receiver(:,1) ).^2+( pos(2)+vel(2)*t-receiver(:,2)).^2 )/v0;
Tx_value=(pos(1)+vel(1)*t-receiver(:,1))./(v0^2*T_value);
Tz_value=(pos(2)+vel(2)*t-receiver(:,2))./(v0^2*T_value);


% a and b  Gauss
a(1)= vel(1)*(sum( (T_value-timeobs+pos(3)).*Tx_value )) / sigma_r^2 + gamma(1);
a(2)= vel(2)*(sum( (T_value-timeobs+pos(3)).*Tz_value )) / sigma_r^2 + gamma(2);
a(3)= vel(3)*(sum( (T_value-timeobs+pos(3)).*1      ))   / sigma_r^2 + gamma(3);
b(1)=vel(1)*(sum((Tx_value*vel(3))))   / sigma_r^2;
b(2)=vel(2)*(sum((Tz_value*vel(3))))   / sigma_r^2;
b(3)=vel(3)*(size(timeobs,1)*vel(3))   / sigma_r^2;
b=abs(b);

M_value=max(0,a+b*t)+gamma;

tau = zeros(dim);
%  a >= 0 and a < 0
    idx_pos = (a >= 0);
    if any(idx_pos)
        ap = a(idx_pos);
        yp = s(idx_pos);
        bp=  b((idx_pos));
        tau(idx_pos) = (-ap + sqrt(ap.^2 + 2*bp.*yp)) ./ bp;
    end
    if any(~idx_pos)
        an = a(~idx_pos);
        yn = s(~idx_pos);
        bn=  b(~idx_pos);
        t0 = -an ./ bn;
        tau(~idx_pos) = t0 + sqrt(2*yn ./ bn);
    end

end