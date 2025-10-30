%calculation of the m
function m_value=rate(t,pos,vel,timeobs,sigma_r,gamma,receiver,v0)
dim=size(pos);
m_value=zeros(dim(1),1);

T_value=sqrt( ( pos(1)+vel(1)*t-receiver(:,1) ).^2+( pos(2)+vel(2)*t-receiver(:,2)).^2 )/v0;
Tx_value=(pos(1)+vel(1)*t-receiver(:,1))./(v0^2*T_value);
Tz_value=(pos(2)+vel(2)*t-receiver(:,2))./(v0^2*T_value);

m_value(1)=sum((T_value+pos(3)+vel(3)*t-timeobs).*Tx_value)/sigma_r^2;
m_value(1)=max([0,vel(1)*m_value(1)])+gamma(1);

m_value(2)=sum((T_value+pos(3)+vel(3)*t-timeobs).*Tz_value)/sigma_r^2;
m_value(2)=max([0,vel(2)*m_value(2)])+gamma(2);

m_value(3)=sum((T_value+pos(3)+vel(3)*t-timeobs).*1       )/sigma_r^2;
m_value(3)=max([0,vel(3)*m_value(3)])+gamma(3);

end