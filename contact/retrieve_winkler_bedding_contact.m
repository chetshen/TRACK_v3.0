a=2e-2;
dx=2e-4;
xq=-a:dx:a;
numTimeStep = 1000;
Fn = zeros(numTimeStep,2);
f = zeros(length(xq),numTimeStep);
for timeStep = 1:numTimeStep
    X_w_t = X_w(1,1)+timeStep*inp.solver.deltat*vx;
    [Fn(timeStep+1,:),f(:,timeStep+1)] = winkler_bedding(X_w(1,1), X_w_t, ...
        dis.r(timeStep+1,:),dis.w(timeStep+1,1),Z.irr, geo, mat_trk,inp,-1);
        disp(num2str(timeStep))
end