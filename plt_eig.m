[NodeNumber]=picknodes([5 25], geo.ND); %define the node numbers to be plotted
[~,~,~,nodeRef]=time_history(NodeNumber,dis,vel,acc, mat_trk.activeDof); %  get the reference index of the nodes to be plotted


Vshape=V(nodeRef,:);
Lia=ismember(geo.ND(:,1),NodeNumber(1:size(Vshape,1)));
shape=zeros(length(geo.ND),1);
shape(Lia)=Vshape(:,277);
plot_modeShape(geo,shape,10);
xlim('Auto');