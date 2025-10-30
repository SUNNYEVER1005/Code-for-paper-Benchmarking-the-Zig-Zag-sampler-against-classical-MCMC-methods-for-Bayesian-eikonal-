function [value1,value2]=gradient_T(T,dx,dz)

[xx,zz]=size(T);

gradTx=zeros(size(T));gradTz=zeros(size(T));

ii=2:xx-1;jj=2:zz-1;

gradTx(ii,jj)=(T(ii+1,jj)-T(ii-1,jj))/(2*dx);
gradTz(ii,jj)=(T(ii,jj+1)-T(ii,jj-1))/(2*dz);

gradTx(1,:)=(T(2,:)-T(1,:))/(2*dx);
gradTx(end,:)=(T(end,:)-T(end-1,:))/(2*dx);
gradTz(:,1)=(T(:,2)-T(:,1))/(2*dz);
gradTz(:,end)=(T(:,end)-T(:,end-1))/(2*dz);  

gradTx(ii,1)=(T(ii+1,1)-T(ii-1,1))/(2*dx);gradTx(ii,end)=(T(ii+1,end)-T(ii-1,end))/(2*dx);
gradTz(1,jj)=(T(1,jj+1)-T(1,jj-1))/(2*dz);gradTz(end,jj)=(T(end,jj+1)-T(end,jj-1))/(2*dz);

value1=gradTx;value2=gradTz;

end