clc;clear;close all;
Num_User=1;
Tx_antBS=8;
RIS_Lnum=8;
Tx_antRIS=32/RIS_Lnum;

power_num=21;
ratemar2=zeros(3,power_num);
flagmar=zeros(2,power_num);
iter_num=3e2;
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
            [ee,rate]=singleuserEEobj(thetamar0,P_max,xonoffini,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
                PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
            
            [ee3,rate0,rate3]=singleuserAF(P_max,xonoffini,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
                PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,P_AFtran);
            
            
            %if isreal(rate)==1
            flagmar (1,p_num)=...
                flagmar (1,p_num)+1;
            ratemar2(1,p_num)=...
                ratemar2(1,p_num)+rate;
            ratemar2(3,p_num)=...
                ratemar2(3,p_num)+rate3;
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
            [ee2,rate2]=singleuserEEobj(thetamar0,P_max,xonoffini,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
                PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
            
            flagmar (2,p_num)=...
                flagmar (2,p_num)+1;
            ratemar2(2,p_num)=...
                ratemar2(2,p_num)+rate2;
            % end
        end
    end
end
%%

if flagmar(1,1)~=0
    ratemar2(1,:)=ratemar2(1,:)./flagmar(1,:);
    ratemar2(3,:)=ratemar2(3,:)./flagmar(1,:);
end

if flagmar(2,1)~=0
    
    ratemar2(2,:)=ratemar2(2,:)./flagmar(2,:);
    
end



%%


figure(3);
xx=powermin:p_delta:powermin+p_delta*power_num-p_delta;
plot(xx, ratemar2(1,:),'-o',xx, ratemar2(2,:),'-s',...
    xx,ratemar2(3,:),'-^','LineWidth',2);
legend('DRIS','CRIS','AFR');
grid on;
xlabel('Maximal transmit power {\it P}_{max} (dBm)','fontsize',12);
ylabel('Sum-rate (Mbps)','fontsize',12);
set(gca,'FontSize',12);


%ylabel('Energy efficiency (Mbits/Joule)','fontsize',12);

