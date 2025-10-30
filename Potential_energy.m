%% U
function value=Potential_energy(theta,timeobs,receiver,v0)  


 % uniform model under posterior distribution (8)


time=sqrt( ( theta(1)-receiver(:,1) ).^2+( theta(2)-receiver(:,2)).^2 )/v0;
sigma_r=0.01;
 
value=(1/sigma_r^2)*0.5*sum(  (time-(timeobs-theta(3))).^2  );

  
end
    
 

    
       
    
        
    
