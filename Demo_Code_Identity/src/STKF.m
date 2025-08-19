function kff=STKF(kf,u,z,nu)

% kf : the kalman fitler instance 
Q=kf.Q;
R=kf.R;
F=kf.F;
H=kf.H;
G=kf.G;
len=kf.len;
n=kf.n;
m=kf.m;
if(nargin==3)
nu = 4;          % Degrees of freedom
end
tau2 = 1;        % Scale parameter (variance)
% stored state
statef_=zeros(n,len);
statef=zeros(n,len);
covf_=zeros(n,n,len);
covf=zeros(n,n,len);
for i=1:len
    if(i==1)
        x=kf.x0;
        P=kf.P0;
    end
    % prediction 
    x_=F*x+G*u(:,i);
    P_=F*P*F'+Q;
    %
    %bp = chol(P_,'lower') ;
    br = chol(R,'lower') ;
    cnt=2;
    num=cnt;
    while(num>0)
        %  
        if(num==cnt)
          x_tlast=x_; 
        else  
          x_tlast=x_t; 
        end
        num=num-1;
        %dp= bp\x_;
        z_=z(:,i);
        dr= br\z_;
        %wp= bp\x_tlast;
        wr= inv(br)*H*x_tlast;
        %ep=dp-wp;
        er=dr-wr;
        %  P_ and R
        %Cy=diag(1/(tau2*(1+er.*er./(nu*tau2))));
        Cy=diag(nu/(tau2*nu+er.*er));
        R_1=br/Cy*br';
        K_1=P_*H'/(H*P_*H'+R_1);
        % obtain the new measurement
        %X_t=X_+K_1*(z(:,i)-H*X_); 
        x_t=x_+K_1*(z_-H*x_);
        xe(cnt-num,i)=norm(x_t-x_tlast)/(norm(x_tlast)+0.00001);
        % stored data for inspectation
        if(num==cnt-1)
        er_matrix(:,i)=er;
        end
        if(num==cnt-2)
        %ep_matrix(:,i)=ep;
        end
        if(xe(cnt-num,i)<0.01)
            break
        end
    end 
    lambda=1/Cy;
    % update
    x=x_t;
    P=(eye(2)-K_1*H)*P_*(eye(2)-K_1*H)'+K_1*R*K_1';
    % store the data
    statef_(:,i)=x_';
    statef(:,i)=x';
    covf_(:,:,i)=P_;
    covf(:,:,i)=P;
    Lambda(i)=lambda;
    iter(i)=cnt-num;
end
kff.statef_=statef_;
kff.statef=statef;
kff.covf_=covf_;
kff.covf=covf;
kff.iter=iter;
kff.er_matrix=er_matrix;
kff.Lambda=Lambda;

end