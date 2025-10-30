function [value1,value2,value3]=traveltime_continous(theta,T_gridip,Tx_gridip,Tz_gridip)  %对输入数组插值



%%  T_m(x_{source},z_{source})
time=zeros(size(T_gridip,1),1);Tx=time;Tz=Tx;  % size(receiver,1)

    for ii=1:size(time,1)
        time(ii)=T_gridip{ii}(theta(1),theta(2));
        Tx(ii)=Tx_gridip{ii}(theta(1),theta(2));
        Tz(ii)=Tz_gridip{ii}(theta(1),theta(2));
    end
    value1= time ;
    value2=Tx;
    value3=Tz;

end

















% %%  T_m(x_{source},z_{source})
% time=zeros(size(T,3),1);Tx=time;Tz=Tx;  % size(receiver,1)
% 
%     for ii=1:size(time,1)
%         time(ii)=matrix_interpolate_function(T(:,:,ii), theta(1),theta(2),x_interval,z_interval,dx,dz);
%         Tx(ii)=matrix_interpolate_function(gradTx(:,:,ii), theta(1),theta(2),x_interval,z_interval,dx,dz);
%         Tz(ii)=matrix_interpolate_function(gradTz(:,:,ii), theta(1),theta(2),x_interval,z_interval,dx,dz);
%     end
%     value1= time ;
%     value2=Tx;
%     value3=Tz;
% 
% end