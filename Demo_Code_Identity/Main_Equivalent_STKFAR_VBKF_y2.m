function Main_Equivalent_STKFAR_VBKF_y2()
% Brief: This code demonstrates the identify of STKF and VBKF_fixed
% Details:
%    None
% 
% Syntax:  
%     Main_Equivalent_STKFAR_VBKF_y2()
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
% Created:                         19-Aug-2025 17:04:17
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
    CAS=2;
    if(CAS==1)
    if(abs(randn(1))<1.96)
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
    amp(i)=(2*sin(2*pi*0.02*t(i)))^2+1;
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
kff=kf_forward(kf,u,z);
% adaptive kalman filter
vbkff=akf_forward_vb(kf,u,z,0.99);
%
kf.nu_p=[10^8 10^8]';
kf.nu_r=[100]';
rhop=1*ones(2,1);
rhor=0.99;
stakf=student_forward_adap(kf,u,z,rhop,rhor);

%%
kff.err=rms(kff.statef-state,2);
vbkff.err= rms(vbkff.statef-state,2);  %
stakf.err= rms(stakf.statef-state,2);  %

fprintf('kf: %.3f, %.3f \r\n', kff.err);
fprintf('vbkf:  %.3f, %.3f \r\n', vbkff.err);
fprintf('rvkf:  %.3f, %.3f \r\n', stakf.err);


% 

figure
subplot(2,1,1)
grid on
plot(t,vbkff.statef(1,:)-state(1,:),'LineWidth',1.2,'Marker','+','MarkerIndices',1:10:length(t))
hold on
plot(t,stakf.statef(1,:)-state(1,:),'LineWidth',1.2)
plot(t,kff.statef(1,:)-state(1,:),'LineWidth',1.2)

legend('VBKF-fixed','STKF','KF','Orientation','horizontal','Interpreter','latex')
xlabel([])
ylabel('$x_1$','Interpreter','latex')
set(gca,'fontsize',16)
subplot(2,1,2)
grid on
plot(t,vbkff.statef(2,:)-state(2,:),'LineWidth',1.2,'Marker','+','MarkerIndices',1:10:length(t))
hold on
plot(t,stakf.statef(2,:)-state(2,:),'LineWidth',1.2)
plot(t,kff.statef(2,:)-state(2,:),'LineWidth',1.2)

grid on;
hold off;
set(gcf,'Position',[100 100 700 600]);
set(gca,'fontsize',16)
xlabel('time (s)','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')

mean(stakf.iter)

%
figure;
% Solid lines
plot(t,stakf.Lambda(3,:)'*0.1, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', 'blue');
hold on;
% Dashed lines
plot(t,vbkff.Lambda'*0.1, 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'red','Marker','square','MarkerIndices',1:200:length(t),'MarkerSize',15,'MarkerEdgeColor','auto');
plot(t,amp.*amp*0.1,'Color',[0 0.7 0.2],'LineWidth',1.5)
% Legend
legend('STKF-AR1 $\sigma^2$', 'VBKF $\sigma^2$', 'Truth',...
       'Interpreter', 'latex','Orientation','horizontal');

% Optional: add labels and title
xlabel('Time (s)','Interpreter', 'latex');
ylabel('Value','Interpreter', 'latex');
grid on;
hold off;
set(gca,'fontsize',16)
set(gcf,'Position',[100 100 700 600]);


end