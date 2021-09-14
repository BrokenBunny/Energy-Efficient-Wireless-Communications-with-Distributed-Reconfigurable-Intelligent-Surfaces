function [xonoff,ee]=singleuserxonoff(thetamar0,power0,xonoff0,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS)

iter=0;
xonoffnew=xonoff0;
ee0=singleuserEEobj(thetamar0,power0,xonoff0,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
 
while iter<=RIS_Lnum
    iter=iter+1;
    
    flag=0;
    for l_RIS=1:RIS_Lnum
        if xonoff0(l_RIS)==1
            xtemp=xonoff0;
            xtemp(l_RIS)=0;
            eenew=singleuserEEobj(thetamar0,power0,xtemp,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
                PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
            if eenew>ee0
                flag=1;
                ee0=eenew;
                xonoffnew=xtemp;
            end
        end
    end
    if flag==1
        xonoff0=xonoffnew;
    else
        break;
    end
end

xonoff=xonoffnew;
ee=ee0;
tes=0;



