 function [dis_out,vel_out,acc_out,index,nodeRef]=time_history(nodeNumber, dis, vel, acc, activeDof, dofID)
if nargin < 6
    dofID=1;
end
dof=2*(nodeNumber-1)+dofID;
[index,loc]=ismember(activeDof,dof);
loc1=loc(loc(:,1)~=0,1);
nodeRef=zeros(length(loc1),1);
for i=1:length(loc1)
nodeRef(i,1)=(dof(loc1(i))-dofID)/2+1;
end

dis_out=dis(:,index);
vel_out=vel(:,index);
acc_out=acc(:,index);
 end
 