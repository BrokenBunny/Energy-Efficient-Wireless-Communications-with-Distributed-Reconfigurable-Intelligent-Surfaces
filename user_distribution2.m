%3-tier heteerogeneous network
%% System parameters
function [Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,...
    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS]=...
    user_distribution2(Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS)
Bandwidth=1;%MHz
Noise=-104;%dBm  -174dBm/Hz
Noise=10^(Noise/10)/10^3;
%% Area: Square
lengA=200;  % m

%power limits
P_max=50; %dBm
P_max=10^(P_max/10)/10^3;
P_k=10;%dBm
P_k=10^(P_k/10)/10^3;
P_R=10;%dBm
P_R=10^(P_R/10)/10^3;
P_A=10;%dBm
P_A=10^(P_A/10)/10^3;
P_B=3;%dBW;
P_B=10^(P_B/10);

mu=1/.8;
%% Users

BS_loc=[lengA/2,lengA/2];

RISloc=zeros(RIS_Lnum,2);
for l_RIS=1:RIS_Lnum
    RISloc(l_RIS,:)=[cos(pi*2*l_RIS/RIS_Lnum),sin(pi*2*l_RIS/RIS_Lnum)]*lengA/3+[lengA/2,lengA/2];
end

User_loc=rand(Num_User,2)*lengA;

dismin=10;%m


%%
alpha=0.2;%-3.53;
beta=3.76;
degrade=1e6;
Distance_UserBS=max(sqrt(sum((repmat(BS_loc,Num_User,1)-User_loc).^2,2)),dismin);
PathLoss_UserBS1=sqrt(10^(alpha)./Distance_UserBS.^(beta));
PathLoss_UserBS=repmat(PathLoss_UserBS1,1,Tx_antBS);
PathLoss_UserBS = PathLoss_UserBS.*1/sqrt(2).*(randn(Num_User,Tx_antBS)+1i*randn(Num_User,Tx_antBS))/degrade;
%%

Distance_UserRIS=zeros(Num_User,RIS_Lnum);
for l_RIS=1:RIS_Lnum
    Distance_UserRIS(:,l_RIS)=max(sqrt(sum((repmat(RISloc(l_RIS,:),Num_User,1)-User_loc).^2,2)),dismin);
end
PathLoss_UserRIS1=sqrt(10^(alpha)./Distance_UserRIS.^(beta));
PathLoss_UserRIS=zeros(Num_User,RIS_Lnum,Tx_antRIS);

for k_user=1:Num_User
    for l_RIS=1:RIS_Lnum
        PathLoss_UserRIS(k_user,l_RIS,:)=repmat(PathLoss_UserRIS1(k_user,l_RIS),Tx_antRIS,1)...
            .*1/sqrt(2).*(randn(Tx_antRIS,1)+1i*randn(Tx_antRIS,1));
    end
end

%%
Distance_RISBS=zeros(RIS_Lnum,1);
for l_RIS=1:RIS_Lnum
    Distance_RISBS(l_RIS,1)=max(sqrt(sum((RISloc(l_RIS,:)-BS_loc).^2)),dismin);
end
PathLoss_RISBS1=sqrt(10^(alpha)./Distance_RISBS.^(beta));
PathLoss_RISBS=zeros(RIS_Lnum,Tx_antRIS,Tx_antBS);


for l_RIS=1:RIS_Lnum
    PathLoss_RISBS(l_RIS,:,:)=repmat(PathLoss_RISBS1(l_RIS,1),Tx_antRIS,Tx_antBS)...
        .*1/sqrt(2).*(randn(Tx_antRIS,Tx_antBS)+1i*randn(Tx_antRIS,Tx_antBS));
end

ttt=0;


%% Date save
%save('SystemData2miso.mat','p_max','PathLoss_User_BS','Noise');

%% Figre plot

% figure(1)
% plot(BS_loc(1,1),BS_loc(1,2),'^','Markersize',10,'Markerfacecolor',[0,0,0],'MarkerEdgecolor','k');
% hold on;
% plot(RISloc(:,1),RISloc(:,2),'o','Markersize',10,'Markerfacecolor',[0,0,0],'MarkerEdgecolor','k');
% hold on;
% plot(User_loc(:,1),User_loc(:,2),'*k');
% xlabel('km');
% ylabel('km');
% legend('BS','RIS','User');
% xlim([0,lengA]);
% ylim([0,lengA]);