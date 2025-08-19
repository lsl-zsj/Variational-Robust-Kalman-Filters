function robust_adaptive_filters_test_y1()
% Brief: This code demonstrate the identiy of VBKF and STKF-AR 
% Details:
%    None
% 
% Syntax:  
%     robust_adaptive_filters_test_y1()
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
% Created:                         19-Aug-2025 17:17:16
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%

%% 

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
[state,z,d,ampw,ampv]=measurement_generation_y1(kfa,uk);
stateaug=[d;state];

%%
vbkff=akf_forward_vb(kfa,uk,z,0.99);


%% adaptive kalman filter
kfa.nu_p=[10^8 10^8 10^8 10^8 10^8]';
kfa.nu_r=[100 100]';
stakf=student_forward_adap(kfa,uk,z);

figure;
% Solid lines
plot(t,stakf.Lambda(6,:)', 'LineStyle', '-', 'LineWidth', 1.5, 'Color', 'blue');
hold on;
plot(t,stakf.Lambda(7,:)', 'LineStyle', '-', 'LineWidth', 1.5, 'Color', 'red');

% Dashed lines
plot(t,vbkff.Lambda(1,:)'*10, 'LineStyle', '-.', 'LineWidth', 1, 'Color', 'black','Marker','square','MarkerIndices',1:200:length(t),'MarkerSize',15,'MarkerEdgeColor','auto');
plot(t,vbkff.Lambda(2,:)'*10, 'LineStyle', '-.', 'LineWidth', 1, 'Color', 'magenta','Marker','+','MarkerIndices',1:200:length(t),'MarkerSize',15,'MarkerEdgeColor','auto');

% Legend
legend('STKF-AR $\sigma_1^2$', 'STKF-AR $\sigma_2^2$', 'VBKF $\sigma_1^2$', 'VBKF $\sigma_2^2$', ...
       'Interpreter', 'latex');

% Optional: add labels and title
xlabel('Time (s)','Interpreter', 'latex');
ylabel('Value','Interpreter', 'latex');
grid on;
hold off;
set(gca,'fontsize',16)

%% error performance
vbkff.err= rms(vbkff.statef-stateaug,2);  %
stakf.err= rms(stakf.statef-stateaug,2);  %

fprintf('vbkf: %.3f, %.3f, %.3f, %.3f, %.3f\n', vbkff.err);
fprintf('rvkf: %.3f, %.3f, %.3f, %.3f, %.3f\n', stakf.err);



end