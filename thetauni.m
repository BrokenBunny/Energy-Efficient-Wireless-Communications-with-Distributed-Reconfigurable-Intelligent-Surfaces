function theta=thetauni(thetavecupd)

theta=thetavecupd;
Num=length(thetavecupd);

for i_num=1:Num
    theta(i_num)=thetavecupd(i_num)/norm(thetavecupd(i_num));    
end

test=0;