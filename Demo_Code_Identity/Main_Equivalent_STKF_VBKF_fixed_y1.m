function Main_Equivalent_STKF_VBKF_fixed_y1()
% Brief: This code demonstrates the identify of STKF and VBKF_fixed
% Details:
%    None
% 
% Syntax:  
%     Main_Equivalent_STKF_VBKF_fixed_y1()
% 
% Inputs:
%    None
% 
% Outputs:
%    None
% 
% Example: 
%    None
% 
% See also: None

% Author:                          ShileiLi
% Email:                           slidk@connect.ust.hk
% Created:                         19-Aug-2025 17:04:50
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%

clear all
addpath(genpath(pwd));
dt=0.01;
R=0.1;
G=[0.5 * dt ^2;
     dt];
Q=G*G';
F=[1 dt;
     0 1];
H=[1 0];
tsim = 50;
len=tsim/dt;
x    = zeros(2, len);
u   = zeros(1, len);            %  input
u_n   = zeros(1, len);          %  input
z    = zeros(1, len);           % measurement at hz
t    = 0:dt:tsim-dt;            % time
x0 = [0;0];
u0=0;
noise=zeros(len,1);
amp=ones(1,len);
indexn=zeros(1,len);
%% state and measurement generation
for i = 1:len
    % process
    u(:,i) = 0;    
    u_n(:,i) = u(:,i)+ randn(1); % this actually is the sequence of 0 to N-1
    if(i==1)
    x(:,i) =F*x0+G*u0;
    else
    x(:,i) = F * x(:,i-1) + G * u_n(:,i);   % Generate truth  u(:,i) actually is u(:,i-1)
    end
    % impulsive noise 
    CAS=1;
    if(CAS==1)
    if(rand<0.95)
        noise(i)=sqrt(R)*randn(1,1);
        indexn(i)=1;
        amp(i)=1;
    else
        noise(i)=10*sqrt(R)*randn(1,1);
        indexn(i)=0;
        amp(i)=10;
    end
    elseif(CAS==2)
    %% slow-varying noise
    amp(i)=abs(2*sin(2*pi*0.05*t(i)))+1;
    % amp(i)=1;
    noise(i)=amp(i)*sqrt(R) *randn(1);
    end
    z(:,i) = H*x(:,i) +  noise(i); 
end

%
state=x;
%
x0=[0;0];
P0 = diag([0.001;0.001]);
kf.Q=Q;
kf.R=R;
kf.F=F;
kf.H=H;
kf.G=G;
kf.x0=x0;
kf.P0=P0;
kf.len=len;
kf.n=2;
kf.m=1;
% kf forward pass
u_n=0*randn(1, len); 

kff=kf_forward(kf,u_n,z);
% adaptive kalman filter
stukf=STKF(kf,u,z);
vbkff=VBKF_fixed(kf,u_n,z);

% 
akf.errpos=rms(vbkff.statef(1,:)-state(1,:));
akf.errvel=rms(vbkff.statef(2,:)-state(2,:));
kf.errpos=rms(kff.statef(1,:)-state(1,:));
kf.errvel=rms(kff.statef(2,:)-state(2,:));
stkf.errpos=rms(stukf.statef(1,:)-state(1,:));
stkf.errvel=rms(stukf.statef(2,:)-state(2,:));

% 
fprintf('vbkf-fixed: %.3f,%.3f \r',akf.errpos,akf.errvel);
fprintf('kf: %.3f ,%.3f \r',kf.errpos,kf.errvel);
fprintf('stkf: %.3f ,%.3f \r',stkf.errpos,stkf.errvel);
% 


figure
subplot(2,1,1)
box on 
grid on
hold on
plot(t,vbkff.statef(1,:)-state(1,:),'LineWidth',1.2,'Marker','+','MarkerIndices',1:10:length(t))
plot(t,stukf.statef(1,:)-state(1,:),'LineWidth',1.2)
plot(t,kff.statef(1,:)-state(1,:),'LineWidth',1.2)
legend('VBKF-fixed','STKF','KF','Orientation','horizontal','Interpreter','latex')
xlabel([])
ylabel('$x_1$','Interpreter','latex')
set(gca,'fontsize',16)

subplot(2,1,2)
box on
grid on
hold on
plot(t,vbkff.statef(2,:)-state(2,:),'LineWidth',1.2,'Marker','+','MarkerIndices',1:10:length(t))
plot(t,stukf.statef(2,:)-state(2,:),'LineWidth',1.2)
plot(t,kff.statef(2,:)-state(2,:),'LineWidth',1.2)
xlabel('time (s)','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'fontsize',16)
set(gcf,'Position',[100 100 700 600]);


mean(stukf.iter)

end