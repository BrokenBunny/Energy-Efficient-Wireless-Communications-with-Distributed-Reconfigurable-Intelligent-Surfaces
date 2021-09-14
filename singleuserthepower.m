function [thetamar,power,ee0,ee1,gainvalue,exitflag]=singleuserthepower(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoff,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini)
exitflag=0;
if Num_User==1
    %% theta optimization
    Qnum=sum(xonoff)*Tx_antRIS;
    thetavec=zeros(Qnum,1);
    g1vec=zeros(Tx_antBS,1);
    g1vec(:,1)=PathLoss_UserBS(1,:);
    userris=zeros(Tx_antRIS,1);
    risbs=zeros(Tx_antRIS,Tx_antBS);
    U1mar=zeros(Qnum,Tx_antBS);
    flag=0;
    for l_RIS=1:RIS_Lnum
        userris(:,1)=PathLoss_UserRIS(1,l_RIS,:);
        risbs(:,:)=PathLoss_RISBS(l_RIS,:,:);%Tx_antRIS,Tx_antBS
        if xonoff(l_RIS)==1
            U1mar(flag*Tx_antRIS+1:(flag+1)*Tx_antRIS,:)=diag(userris')*risbs;%;Tx_antRIS,Tx_antRIS;Tx_antRIS,Tx_antBS
            flag=flag+1;
        end
    end
    
    flag=0;
    thetavecini=zeros(Qnum,1);
    for l_RIS=1:RIS_Lnum
        if xonoff(l_RIS)==1
            thetavecini(flag*Tx_antRIS+1:(flag+1)*Tx_antRIS,1)=thetamarini((l_RIS-1)*Tx_antRIS+1:(l_RIS)*Tx_antRIS,1);
            flag=flag+1;
        end
    end
    itermax=2e1;
    epsilon=1e-5;
    gainvalue=zeros(itermax+1,1);
    gainvalue(1)=norm(g1vec+U1mar'*thetavecini);
    %%
    for iter=1:itermax
        thetavecupd=U1mar*(g1vec+U1mar'*thetavecini);
        thetavecupd=thetauni(thetavecupd);
        gainvalue(iter+1)=norm(g1vec+U1mar'*thetavecupd);
        if abs(gainvalue(iter+1)-gainvalue(iter))/gainvalue(iter)<epsilon&&iter>10
            thetavec=conj(thetavecupd);
            break;
        end
        thetavecini=thetavecupd;
    end
    %%
    thetavec=conj(thetavecupd);
    
    %     %%
    %     [Vec,Dnum] = eig(U1mar*U1mar');% Vec(:,qnum)
    %     %Dnum=diag(Dnum);
    %     pathsquare=0;
    %     for qnum=1:Qnum
    %         thetavectemp=Vec(:,qnum);
    %         temp=norm(g1vec'+thetavectemp'*U1mar);
    %         if temp>pathsquare
    %             pathsquare=temp;
    %             thetavec=conj(thetavectemp);
    %         end
    %     end
    
    thetamar=zeros(RIS_Lnum*Tx_antRIS,1);
    flag=0;
    for l_RIS=1:RIS_Lnum
        if xonoff(l_RIS)==1
            thetamar((l_RIS-1)*Tx_antRIS+1:(l_RIS)*Tx_antRIS,1)=thetavec(flag*Tx_antRIS+1:(flag+1)*Tx_antRIS,1);
            flag=flag+1;
        end
    end
    %% power optimization
    barg1=norm(g1vec'+thetavec.'*U1mar)^2/Noise;
    P_min=(2^(Ratemin/Bandwidth)-1)/barg1;
    if P_min>P_max
        exitflag=1;
    end
    
    P_0=P_k+P_R*Qnum+P_B;
    
    %power=(-mu+P_0*barg1)/(mu*barg1*lambertw(1/exp(1)*(-mu+P_0*barg1)/mu))-1/barg1;
    iter_max=2e1;
    power=P_max;
    for iter=1:iter_max
        lambda=Bandwidth*log2(1+barg1*power)/(mu*power+P_0);
        power=Bandwidth/log(2)/lambda/mu-1/barg1;
        power=min(max(power,P_min),P_max);
    end
    
    te=barg1*(mu*power+P_0)-mu*(1+barg1*power)*log(1+barg1*power);
    power=min(max(power,P_min),P_max);
    ee0=Bandwidth*log2(1+barg1*powerini)/(mu*powerini+P_0);
    ee1=Bandwidth*log2(1+barg1*power)/(mu*power+P_0);
    
    %     %%
    %     g1vec_2=PathLoss_UserBS(1,:);
    %     g1vec_2=g1vec_2';
    %     for l_RIS=1:RIS_Lnum
    %         userris(:,1)=PathLoss_UserRIS(1,l_RIS,:);
    %         risbs(:,:)=PathLoss_RISBS(l_RIS,:,:);%Tx_antRIS,Tx_antBS
    %
    %         thetal=thetamar((l_RIS-1)*Tx_antRIS+1:(l_RIS)*Tx_antRIS,1);
    %
    %         g1vec_2=g1vec_2+userris'*(xonoff(l_RIS)*diag(thetal)+(1-xonoff(l_RIS))*eye(length(thetal)))*risbs;
    %     end
    %
    %
    %
    %     barg1_2=norm(g1vec_2)^2/Noise;
    %     P_0_2=P_k+P_R*sum(xonoff)*Tx_antRIS+P_B;
    %     ee=Bandwidth*log2(1+barg1_2*power)/(mu*power+P_0_2);
    %%
else
    thetamar=0;
    power=0;
    ee0=0;
    ee1=0;
end
te=0;