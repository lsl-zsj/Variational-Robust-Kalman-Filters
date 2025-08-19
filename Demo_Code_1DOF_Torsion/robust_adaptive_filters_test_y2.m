function robust_adaptive_filters_test_y2()
% Brief: This code demonstrates the ability of STKF-AR to adapt to
% time-varying process noise.
% Details:
%    None
% 
% Syntax:  
%     robust_adaptive_filters_test_y2()
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
% Created:                         19-Aug-2025 17:15:08
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%

clear all
addpath(genpath(pwd));  % 将当前目录加入路径
load('sys.mat')

%
dt=0.01;
F=sys.Ad;
G1=sys.Bd1;
G2=sys.Bd2;
H=sys.Cd;

% kd-dob 
Fa=[1 zeros(1,4);       % state transfer matrix
    G2,F];
Ga=[0;G1];              % input matrix
Ha=[zeros(2,1),H];      % observation matrix 

Vwd=0.01;
Vwn=0.01;
Vvn=0.1;


Q=[1;1;(1/dt)^2;(1/dt)^2];
% Q=[1;1;1;1];
Qa=diag([Vwd;Vwn.*Q]);

Ra=Vvn*eye(2);              % measurement covariance
x0a=[0;0;0;0;0];        % initial guess
P0a=eye(5);             % initial covariance

tsim = 20;
len=tsim/dt;            % step number
t    = 0:dt:tsim-dt;          % time

kfa.dt=dt;
kfa.F=Fa; 
kfa.H=Ha;
kfa.Ga=Ga; 
kfa.Q=Qa;
kfa.R=Ra;
kfa.x0=x0a;
kfa.P0=0.01*P0a;
kfa.len=len;
kfa.n=5;
kfa.m=2;

%%
% ground truth and measurement generation
uk=0*ones(1,len);
[state,z,d]=measurement_generation_y2(kfa,uk);
stateaug=[d;state];

%%
vbkff=akf_forward_vb(kfa,uk,z,0.98);


%% adaptive kalman filter
kfa.nu_p=[100 100 100 100 100]';
kfa.nu_r=[10^8 10^8]';
rho_p=[0.98;0.98;0.98;0.98;0.98];
rho_r=1*ones(2,1);
stakf=student_forward_adap(kfa,uk,z,rho_p,rho_r);
 

%% KF-DOB
kfa.Q=Qa;
kfdob=kf_dob_forward(kfa,uk,z);

%% RBKF1: correntropy (0,infty)
kfa.sigma_p=[2 2 2 2 2]';
kfa.sigma_r=[10^8 10^8]';
rbkf1=rbkf1_forward(kfa,uk,z);
%% RBKF2:  (0,2)
kfa.sigma_p=[0.5 0.5 0.5 0.5 0.5]';
kfa.sigma_r=[1.999 1.999]';
rbkf2=rbkf2_forward(kfa,uk,z);

%% RBKF3: (0,infty)
kfa.sigma_p=[1 1 1 1 1]';
kfa.sigma_r=[10^8 10^8]';
rbkf3=rbkf3_forward(kfa,uk,z);




% figure
% plot(t,stakf.Lambda(2:3,:)','LineWidth',1)
% xlabel('Time (s)','Interpreter', 'latex');
% ylabel('Value','Interpreter', 'latex');
% grid on;
% hold off;
% set(gca,'fontsize',16)
% legend('VRKF ${\sigma_{\theta_1}}^2$', 'VRKF ${\sigma_{\theta_2}}^2$','Interpreter', 'latex');

% hold on
% plot(vbkff.Lambda'/Ra)
% legend

%% error performance
vbkff.err= rms(vbkff.statef-stateaug,2);  %kf-dob
stakf.err= rms(stakf.statef-stateaug,2);  %rbkf1
kfdob.err= rms(kfdob.statef-stateaug,2);  %rbkf1
rbkf1.err= rms(rbkf1.statef-stateaug,2);  %rbkf1
rbkf2.err= rms(rbkf2.statef-stateaug,2);  %rbkf1
rbkf3.err= rms(rbkf3.statef-stateaug,2);  %rbkf1

%% visualization

figure
box on
hold on
plot(t,vbkff.statef(2:3,:)-stateaug(2:3,:),'LineWidth',1)
plot(t,stakf.statef(2:3,:)-stateaug(2:3,:),'LineWidth',1)
xlabel('Time (s)','Interpreter', 'latex');
ylabel('Angle error','Interpreter', 'latex');
grid on;
hold off;
set(gca,'fontsize',16)
legend('VBKF ${{\theta_1}}^2$', 'VBKF ${{\theta_2}}^2$',...
    'STKF-AR ${{\theta_1}}^2$', 'STKF-AR ${{\theta_2}}^2$','Interpreter', 'latex');


%% error performance
fprintf('vbkf: &%.3f  &%.3f  &%.3f  &%.3f  &%.3f\n', vbkff.err);
fprintf('stkf-ar: &%.3f  &%.3f  &%.3f  &%.3f  &%.3f\n', stakf.err);
fprintf('kfdob: &%.3f  &%.3f  &%.3f  &%.3f  &%.3f\n', kfdob.err);
fprintf('rbkf1: &%.3f  &%.3f  &%.3f  &%.3f  &%.3f\n', rbkf1.err);
fprintf('rbkf2: &%.3f  &%.3f  &%.3f  &%.3f  &%.3f\n', rbkf2.err);
fprintf('rbkf3: &%.3f  &%.3f  &%.3f  &%.3f  &%.3f\n', rbkf3.err);



end