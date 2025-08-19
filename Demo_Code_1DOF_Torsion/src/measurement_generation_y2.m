function [state,zz,dd]=measurement_generation_y2(kfa,uk)


% obtain the system dynamcis


len=kfa.len;
F=kfa.F(2:5,2:5);                % state transfer matrix
G2=kfa.F(2:5,1);                 % for disturbance
G1=kfa.Ga(2:5);                  % for input
H=kfa.H(:,2:5);                % observation matrix
Q=kfa.Q(2:5,2:5);
Qd=kfa.Q(1,1);
R=kfa.R;

%rng(42)
dk=sqrt(Qd)*randn(1,len);
wk=sqrt(Q)*randn(4,len);           % for process noise
vk1=sqrt(R)*randn(2,len);
vk2=sqrt(R)*randn(2,len);
x    = zeros(4, len);
d   = zeros(1, len);          % control input
z    = zeros(2, len);         % measurement at hz
x0 = zeros(4,1);

ampw=ones(len,1);
ampv=ones(len,1);

for i = 1:len
    % process disturbance generation
    if(i>800&&i<1200)
        d(:,i) = 0;
        dk(i)=sqrt(100*Qd)*randn(1);
        wk(:,i)=sqrt(100*Q)*randn(4,1);
    else
        d(:,i) = 0;   
        dk(i)=sqrt(Qd)*randn(1);
        wk(:,i)=sqrt(Q)*randn(4,1);
    end


    d(:,i)= d(:,i) + dk(i);
    % real state generation
    if(i==1)
    x(:,i) =F*x0 + G1*uk(:,i) + G2*d(:,i)+ wk(:,i); % no input and disturbance
    else
    x(:,i) = F * x(:,i-1) + G1*uk(:,i)+ G2*d(:,i) + wk(:,i);   % Generate truth  u(:,i) actually is u(:,i-1)
    end
    % measurement
    vk1(:,i)=sqrt(R)*randn(2,1);
    vk2(:,i)=sqrt(R)*randn(2,1);
    if(i<=5000)
    z(:,i) = H*x(:,i) +   vk1(:,i);  
    else
    z(:,i) = H*x(:,i) +   vk2(:,i); 
    end


end
state=x;
zz=z;
dd=d;




end