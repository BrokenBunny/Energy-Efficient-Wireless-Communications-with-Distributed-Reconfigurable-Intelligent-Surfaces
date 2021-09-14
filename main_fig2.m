clc;clear;close all;
Num_User=1;
Tx_antBS=8;
RIS_Lnum=8;
Tx_antRIS=32/RIS_Lnum;

power_num=21;
eemar2=zeros(5,power_num);
flagmar=zeros(2,power_num);
iter_num=5e1;
temp=zeros(2,iter_num);
powermin=10;
p_delta=2;
P_AFtran=0;
P_AFtran=10^(P_AFtran/10);
%%
for iter=1:iter_num
    [Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
        PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS]=...
        user_distribution2(Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
    P_max=powermin;
    P_max=10^(P_max/10)/1e3;
    Ratemin=ones(Num_User,1);% 1Mbits
    powerini=P_max;
    xonoffini=ones(RIS_Lnum,1);
    thetamarini=ones(RIS_Lnum*Tx_antRIS,1);
    [thetamar0,power0,ee0,ee1,gainvalue,exitflag1]=singleuserthepower(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
        PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
    
    
    if exitflag1==0
        for p_num=1:power_num
            P_max=powermin+p_num*p_delta-p_delta;
            P_max=10^(P_max/10)/1e3;
            [thetamar,power,xonoff,ee,rate]=singleuseroptmi(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
                PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
            
            [ee3,rate3]=singleuserAF(power,xonoff,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
                PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,P_AFtran);
            
            ee5=singleuseroptmiexh(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
                PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
            %if isreal(rate)==1
            flagmar (1,p_num)=...
                flagmar (1,p_num)+1;
            eemar2(1,p_num)=...
                eemar2(1,p_num)+ee;
            eemar2(3,p_num)=...
                eemar2(3,p_num)+ee3;
            eemar2(5,p_num)=...
                eemar2(5,p_num)+ee5;
            % end
        end
    end
end
%%
RIS_Lnum=1;
Tx_antRIS=32/RIS_Lnum;
for iter=1:iter_num
    [Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
        PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS]=...
        user_distribution2(Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
    P_max=powermin;
    P_max=10^(P_max/10)/1e3;
    Ratemin=ones(Num_User,1);% 1Mbits
    powerini=P_max;
    xonoffini=ones(RIS_Lnum,1);
    thetamarini=ones(RIS_Lnum*Tx_antRIS,1);
    [thetamar0,power0,ee0,ee1,gainvalue,exitflag1]=singleuserthepower(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
        PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
    
    
    if exitflag1==0
        for p_num=1:power_num
            P_max=powermin+p_num*p_delta-p_delta;
            P_max=10^(P_max/10)/1e3;
            [thetamar2,power2,xonoff2,ee2,rate2]=singleuseroptmi(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
                PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
            
            flagmar (2,p_num)=...
                flagmar (2,p_num)+1;
            eemar2(2,p_num)=...
                eemar2(2,p_num)+ee2;
            % end
        end
    end
end
%%

if flagmar(1,1)~=0
    eemar2(1,:)=eemar2(1,:)./flagmar(1,:);
    eemar2(3,:)=eemar2(3,:)./flagmar(1,:);
    eemar2(5,:)=eemar2(5,:)./flagmar(1,:);
end

if flagmar(2,1)~=0
    
    eemar2(2,:)=eemar2(2,:)./flagmar(2,:);
    
end



%%


figure(2);
xx=powermin:p_delta:powermin+p_delta*power_num-p_delta;
plot(xx, eemar2(1,:),'-o',xx, eemar2(2,:),'-s',...
    xx,eemar2(3,:),'-^',xx,eemar2(5,:),'-+','LineWidth',2);
legend('DRIS','CRIS','Relay','EXH-DRIS');
grid on;
xlabel('Maximal transmission power {\it P}_{max} (dBm)','fontsize',12);
ylabel('Average EE (Mbits/Joule)','fontsize',12);
set(gca,'FontSize',12);

