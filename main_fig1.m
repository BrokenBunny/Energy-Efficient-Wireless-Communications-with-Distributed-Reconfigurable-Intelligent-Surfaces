clc;clear;close all;
Num_User=1;
Tx_antBS=8;
RIS_Lnum=8;
Tx_antRIS=32/RIS_Lnum;

[Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS]=...
    user_distribution2(Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS);
Ratemin=ones(Num_User,1);% 1Mbits
powerini=P_max;
xonoffini=ones(RIS_Lnum,1);
thetamarini=ones(RIS_Lnum*Tx_antRIS,1);

[thetamar0,power0,ee0,ee1,objmar1]=singleuserthepower(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
objmar1=objmar1.^2/Noise;
 
thetamarini=randn(RIS_Lnum*Tx_antRIS,1)+1i*randn(RIS_Lnum*Tx_antRIS,1); 
thetamarini=thetauni(thetamarini);
[thetamar0,power0,ee0,ee1,objmar2]=singleuserthepower(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
objmar2=objmar2.^2/Noise;
 thetamarini=randn(RIS_Lnum*Tx_antRIS,1)+1i*randn(RIS_Lnum*Tx_antRIS,1); 
thetamarini=thetauni(thetamarini);
[thetamar0,power0,ee0,ee1,objmar3]=singleuserthepower(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoffini,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,Ratemin,powerini,thetamarini);
objmar3=objmar3.^2/Noise;

x_mar=1:10;
plot(x_mar,objmar1(x_mar),'-o',...
    x_mar,objmar2(x_mar),'-s',...
    x_mar,objmar3(x_mar),'-v',...
    'linewidth',2);
xlim([min(x_mar),max(x_mar)]);
xlabel('Iteration number');
ylabel('Normalized channle gain value');
legend('Random initial solution 1','Random initial solution 2','Random initial solution 3')
set(gca,'fontsize',12);
grid on;