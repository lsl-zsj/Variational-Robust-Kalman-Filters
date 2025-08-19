function [state,zz,dd,amp]=measurement_generation_y4(kfa,uk)

% coexist adaptive and outlier noise in measurement 

% obtain the system dynamcis

len=kfa.len;
F=kfa.F(2:5,2:5);                % state transfer matrix
G2=kfa.F(2:5,1);                 % for disturbance
G1=kfa.Ga(2:5);                  % for input
H=kfa.H(:,2:5);                % observation matrix
Q=kfa.Q(2:5,2:5);
Qd=kfa.Q(1,1);
R=kfa.R;
dt=kfa.dt;
rng(44)
dk=sqrt(Qd)*randn(1,len);
wk=sqrt(Q)*randn(4,len);           % for process noise
vk1=sqrt(R)*randn(2,len);
vk2=sqrt(R)*randn(2,len);
x    = zeros(4, len);
d   = zeros(1, len);          % control input
z    = zeros(2, len);         % measurement at hz
x0 = zeros(4,1);
amp=zeros(len,1);


for i = 1:len
    % process disturbance generation
    if(rand<0.99)
        d(:,i) = 0;
        dk(i)=sqrt(1*Qd)*randn(1);
        wk(:,i)=sqrt(1*Q)*randn(4,1);
    else
        d(:,i) = 0;   
        dk(i)=sqrt(900*Qd)*randn(1);
        wk(:,i)=sqrt(900*Q)*randn(4,1);
    end


    d(:,i)= d(:,i) + dk(i);
    % real state generation
    if(i==1)
    x(:,i) =F*x0 + G1*uk(:,i) + G2*d(:,i)+ wk(:,i); % no input and disturbance
    else
    x(:,i) = F * x(:,i-1) + G1*uk(:,i)+ G2*d(:,i) + wk(:,i);   % Generate truth  u(:,i) actually is u(:,i-1)
    end
    % measurement
    t=i*dt;
    vk1(:,i)=(2*abs(sin(2*pi*0.05*t))^2+1)*sqrt(R)*randn(2,1);
    amp(i)=(2*abs(sin(2*pi*0.05*t))^2+1);
    z(:,i) = H*x(:,i) +   vk1(:,i); 
    
    

end
state=x;
zz=z;
dd=d;



end