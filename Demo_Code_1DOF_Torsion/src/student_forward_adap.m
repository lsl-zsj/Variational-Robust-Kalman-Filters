function kff=student_forward_adap(kf,u,z,rhop,rhor)

% kf : the kalman fitler instance 
Q=kf.Q;
R=kf.R;
F=kf.F;
H=kf.H;
G=kf.Ga;
len=kf.len;
n=kf.n;
m=kf.m;
%
nu_p=kf.nu_p;
nu_r=kf.nu_r;
%
tau2p = ones(n,1);
tau2r = ones(m,1);        % Scale parameter (variance)
% stored state
statef_=zeros(n,len);
statef=zeros(n,len);
covf_=zeros(n,n,len);
covf=zeros(n,n,len);
if(nargin==3)
   rhop=1*ones(n,1);
   rhor=0.99*ones(m,1);
else
end
lambdap=ones(n,1);
lambdar=ones(m,1);
for i=1:len
    if(i==1)
        x=kf.x0;
        P=kf.P0;
    end
    if(i==2001)
        flag=1;
    end
    % prediction 
    x_=F*x+G*u(:,i);
    P_=F*P*F'+Q;
    %
    nu_p_=rhop.*nu_p;
    nu_p=nu_p_+1;
    nu_r_=rhor.*nu_r;
    nu_r=nu_r_+1;
    tau2p=rhop.*tau2p; % priori
    tau2r=rhor.*tau2r; % priori
    bp = chol(P_,'lower') ;
    br = chol(R,'lower') ;
    cnt=4;
    num=cnt;
    while(num>0)
        %  
        if(num==cnt)
          x_tlast=x_; 
        else  
          x_tlast=x_t; 
        end
        num=num-1;
        dp= bp\x_;
        z_=z(:,i);
        dr= br\z_;
        wp= bp\x_tlast;
        wr= br\H*x_tlast;
        ep=dp-wp;
        er=dr-wr;
        %  P_ and R
        Cx=diag(nu_p./(tau2p.*nu_p+ep.*ep)); %
        Cy=diag(nu_r./(tau2r.*nu_r+er.*er)); %
        P_1=bp/Cx*bp';
        R_1=br/Cy*br';
        K_1=P_1*H'/(H*P_1*H'+R_1);
        % obtain the new measurement
        %X_t=X_+K_1*(z(:,i)-H*X_); 
        x_t=x_+K_1*(z_-H*x_);
        xe(cnt-num,i)=norm(x_t-x_tlast)/(norm(x_tlast)+0.00001);

        if(xe(cnt-num,i)<0.000001)
            break
        end
    end 

    x=x_t;
    P=(eye(n)-K_1*H)*P_*(eye(n)-K_1*H)'+K_1*R*K_1';


    W=[inv(bp),zeros(n,m);zeros(m,n),inv(br)]*[eye(n);H];
    vector=diag(W*P*W');
    vp=vector(1:n);
    vr=vector(n+1:end);
    %% update
    if(Cx(1,1)<0.9)
        flagp=1;
    else
        flgar=1;
    end
    
    for ii=1:n
        if(rhop(ii)==1)
        tau2p(ii)=lambdap(ii);
        else
        lambdap(ii)=rhop(ii)*lambdap(ii)+(ep(ii))*(ep(ii)')./nu_p(ii)+vp(ii)./nu_p(ii);
        tau2p(ii)=lambdap(ii);
        end
    end
    for ii=1:m
        if(rhor(ii)==1)
        tau2r(ii)=lambdar(ii);
        else
        lambdar(ii)=rhor(ii)*lambdar(ii)+(er(ii))*(er(ii)')./nu_r(ii)+vr(ii)./nu_r(ii);
        tau2r(ii)=lambdar(ii);
        end
    end
    
    % store the data
    statef_(:,i)=x_';
    statef(:,i)=x';
    covf_(:,:,i)=P_;
    covf(:,:,i)=P;
    LAMbda(:,i)=[lambdap;lambdar];
    NU(:,i)=[nu_p;nu_r];
    TAU(:,i)=[tau2p;tau2r];
end
kff.statef_=statef_;
kff.statef=statef;
kff.covf_=covf_;
kff.covf=covf;
kff.Lambda=LAMbda;
kff.NU=NU;
kff.TAU=TAU;
end