function kff=VBKF_fixed(kf,u,z)
% kf : the kalman fitler instance 
Q=kf.Q;
R=kf.R;
F=kf.F;
H=kf.H;
G=kf.G;
len=kf.len;
n=kf.n;
m=kf.m;
LambdaS=zeros(1,len);
stateS=zeros(n,len);
covS=zeros(n,n,len);
ALPHA=zeros(1,len);
BETA=zeros(1,len);
lambdat=1;
alpha = 2;    % Shape parameter (α)
beta = 2;     % Scale parameter (β) \beta=\alpha+1
N=4;
m=1;
for i=1:len
    if(i==1)
    xx=kf.x0;
    PP=kf.P0;
    end
    %%
    alpha=2;
    beta_=2;
    xx_=F*xx;
    PP_=F*PP*F'+Q;
    br = chol(R,'lower') ;
    for j=1:N
        lambdat=beta./(alpha);
        % R_t=br*diag(lambdat)*br';
        R_t=br*diag(lambdat)*br';
        % update
        K_1=PP_*H'/(H*PP_*H'+R_t);
        xx=xx_+K_1*(z(:,i)-H*xx_);
        PP=(eye(n)-K_1*H)*PP_*(eye(n)-K_1*H)'+K_1*R*K_1';
        for jj=1:m
            err=inv(br)*(z(:,i)-H*xx);
            prr=diag(inv(br)*H*PP*H'*(inv(br))');
            beta(jj) = beta_(jj)+1/2*(err(jj))'*(err(jj))+1/2*prr(jj);
        end
    end
    %
    % store the data
    ALPHA(i)=alpha;
    BETA(i)=beta;
    LambdaS(i)=lambdat;
    stateS(:,i)=xx;
    covS(:,:,i)=PP;
end

kff.Lambda=LambdaS;
kff.statef=stateS;
kff.covf=covS;
kff.ALPHA=ALPHA;
kff.BETA=BETA;

end