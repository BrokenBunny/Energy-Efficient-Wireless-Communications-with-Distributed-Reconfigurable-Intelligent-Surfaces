function ee=singleuseroptmiexh(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini)
itermax=1e1;
eemar=zeros(1,itermax+1);
[thetamar,power,xonoff,ee,rate]=singleuseroptmi(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
eemar(itermax+1)=ee;
for iter=1:itermax
    xonoffini=randn(RIS_Lnum,1);
    xonoffini(xonoffini>=-.5)=1;
    xonoffini(xonoffini<0)=0;
    thetamarini=randn(RIS_Lnum*Tx_antRIS,1)+1i*randn(RIS_Lnum*Tx_antRIS,1);
    thetamarini=thetauni(thetamarini);
    [thetamar,power,xonoff,ee,rate]=singleuseroptmi(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
        PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
    eemar(iter)=ee;
    
end

ee=max(eemar);

