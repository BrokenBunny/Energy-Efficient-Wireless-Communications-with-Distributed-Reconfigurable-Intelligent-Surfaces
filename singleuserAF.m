function [ee,rate,ratemax]=singleuserAF(power2,xonoff,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,P_AFtran)

thetamar=ones(RIS_Lnum*Tx_antRIS,1);

g1vec=zeros(Tx_antBS,1);
g1vec(:,1)=PathLoss_UserBS(1,:);
g1vec=g1vec';
userris=zeros(Tx_antRIS,1);
risbs=zeros(Tx_antRIS,Tx_antBS);
Noisegain=0;
for l_RIS=1:RIS_Lnum
    userris(:,1)=PathLoss_UserRIS(1,l_RIS,:);
    risbs(:,:)=PathLoss_RISBS(l_RIS,:,:);%Tx_antRIS,Tx_antBS
    
    thetal=thetamar((l_RIS-1)*Tx_antRIS+1:(l_RIS)*Tx_antRIS,1);
    temp=norm(risbs);
    temp2=(P_AFtran/norm(userris'));
    g1vec=g1vec+userris'*temp2*risbs;
    Noisegain=Noisegain+temp;
end



barg1=norm(g1vec)^2/Noise/(1+temp*temp2);
P_0=P_k+P_A*sum(xonoff)*Tx_antRIS+P_B+P_AFtran*RIS_Lnum;

P_min=(2^(Ratemin/Bandwidth)-1)/barg1;
iter_max=2e1;
power=P_max;
for iter=1:iter_max
    lambda=Bandwidth*log2(1+barg1*power)/(mu*power+P_0);
    power=Bandwidth/log(2)/lambda/mu-1/barg1;
    power=min(max(power,P_min),P_max);
end

ee=Bandwidth*log2(1+barg1*power)/(mu*power+P_0);
rate=Bandwidth*log2(1+barg1*power);
ratemax=Bandwidth*log2(1+barg1*power2);