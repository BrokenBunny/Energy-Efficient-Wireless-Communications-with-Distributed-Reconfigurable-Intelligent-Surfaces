function [ee,rate]=singleuserEEobj(thetamar,power,xonoff,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS)

g1vec=zeros(Tx_antBS,1);
g1vec(:,1)=PathLoss_UserBS(1,:);
g1vec=g1vec';
userris=zeros(Tx_antRIS,1);
risbs=zeros(Tx_antRIS,Tx_antBS);

for l_RIS=1:RIS_Lnum
    userris(:,1)=PathLoss_UserRIS(1,l_RIS,:);
    risbs(:,:)=PathLoss_RISBS(l_RIS,:,:);%Tx_antRIS,Tx_antBS
    
    thetal=thetamar((l_RIS-1)*Tx_antRIS+1:(l_RIS)*Tx_antRIS,1);
    
    g1vec=g1vec+userris'*(xonoff(l_RIS)*diag(thetal))*risbs;
end



barg1=norm(g1vec)^2/Noise;
P_0=P_k+P_R*sum(xonoff)*Tx_antRIS+P_B;
ee=Bandwidth*log2(1+barg1*power)/(mu*power+P_0);
rate=Bandwidth*log2(1+barg1*power);
