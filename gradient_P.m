%% gradient U

function value=gradient_P(theta,timeobs,receiver,v0)  


% uniform model under posterior distribution (8)

time=sqrt( ( theta(1)-receiver(:,1) ).^2+( theta(2)-receiver(:,2)).^2 )/v0;
sigma_r=0.01;


value=(1/sigma_r^2)*[sum(  ( time-(timeobs-theta(3)) ).*(theta(1)-receiver(:,1))./(v0^2*time) );
           sum(  ( time-(timeobs-theta(3)) ).*(theta(2)-receiver(:,2))./(v0^2*time) );
           sum(    time-(timeobs-theta(3)) )                                        ];



end