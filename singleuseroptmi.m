function [thetamar,power,xonoff,eevalue,rate]=singleuseroptmi(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoff0,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini)

iter=1;
iter_max=1e2;
epsilon=1e-5;
eemar=zeros(1,iter_max*3+1);
thetamar=ones(RIS_Lnum*Tx_antRIS,1);
eemar(1)=singleuserEEobj(thetamarini,powerini,xonoff0,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
powerini=P_max;
while iter<iter_max
    
    [thetamar,power,eemar(iter*3-1),eemar(iter*3)]=singleuserthepower(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
        PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoff0,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
    
    temp=singleuserEEobj(thetamar,power,xonoff0,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
        PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
    
    [xonoff,eemar(iter*3+1)]=singleuserxonoff(thetamar,power,xonoff0,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
        PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
    xonoff0=xonoff;
    powerini=power;
    thetamarini=thetamar;
    if abs(eemar(iter*3+1)-eemar(iter*3))/eemar(iter*3+1)<epsilon
        break;
    end
    iter=iter+1;
end
eevalue=eemar(iter*3);
[temp,rate]=singleuserEEobj(thetamar,power,xonoff,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
        PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
te=0;