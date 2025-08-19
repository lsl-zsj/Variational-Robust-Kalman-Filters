function [state,zz,dd,ampw,ampv]=measurement_generation_y1(kfa,uk)

% obtain the system dynamcis

len=kfa.len;
F=kfa.F(2:5,2:5);                % state transfer matrix
G2=kfa.F(2:5,1);                 % for disturbance
G1=kfa.Ga(2:5);                  % for input
H=kfa.H(:,2:5);                % observation matrix
Q=kfa.Q(2:5,2:5);
Qd=kfa.Q(1,1);
R=kfa.R;

rng(42)
dk=1*sqrt(Qd)*randn(1,len);
wk=sqrt(Q)*randn(4,len);           % for process noise

vk1=sqrt(R)*randn(2,len);
vk2=sqrt(10*R)*randn(2,len);

x    = zeros(4, len);
d   = zeros(1, len);          % control input
z    = zeros(2, len);         % measurement at hz
x0 = zeros(4,1);
amp=0; %

ampw=ones(len,1);
ampv=ones(len,1);

for i = 1:len
    % process disturbance generation
    if(i>200&&i<400)
        d(:,i) = amp;   
    else
        d(:,i) = 0;    
    end
    if(i==201||i==401)
        ampw(i)=10000;
    else
        ampw(i)=1;
    end

    d(:,i)= d(:,i) + dk(i);
    % real state generation
    if(i==1)
    x(:,i) =F*x0 + G1*uk(:,i) + G2*d(:,i)+ wk(:,i); % no input and disturbance
    else
    x(:,i) = F * x(:,i-1) + G1*uk(:,i)+ G2*d(:,i) + wk(:,i);   % Generate truth  u(:,i) actually is u(:,i-1)
    end
    % measurement
    if(i<=500)
    z(:,i) = H*x(:,i) +   vk1(:,i);  
    ampv(i)=1;
    else
    z(:,i) = H*x(:,i) +   vk2(:,i); 
    ampv(i)=1;
    end


end
state=x;
zz=z;
dd=d;




end