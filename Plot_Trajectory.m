
states = in_profile;
out_states = out_profile;
%==========================================================================
%convert curvilinear position to Cartesian
K=size(states,1);
rN = zeros(K,1);
rE = zeros(K,1);
for k=1:K
    [R_N,~] = Radii_of_curvature(states(k,2));
    rN (k,1) = ((states(k,2) ) * (R_N + states(k,4)) * cos(states(k,2)));
    rE (k,1) = ((states(k,3) ) * (R_N + states(k,4)));
end    
%--------------------------------------------------------------------------
%convert curvilinear position to Cartesian
K=size(out_states,1);
rN2 = zeros(K,1);
rE2 = zeros(K,1);
for k=1:K
    [R_N2,~] = Radii_of_curvature(out_states(k,2));
    rN2 (k,1) = ((out_states(k,2) ) * (R_N2 + out_states(k,4)) * cos(out_states(k,2)));
    rE2 (k,1) = ((out_states(k,3) ) * (R_N2 + out_states(k,4)));
end    
%--------------------------------------------------------------------------
%plot 3-D true motion path
fig = figure;
set(fig,'units','normalized');
set(fig,'OuterPosition',[0.05,0.4,0.45,0.6]);
plot3(rN,rE,states(:,4));
hold on;
plot3(rN2,rE2,out_states(:,4));
grid
title('3-D motion path');
xlabel('North (m)');
ylabel('East (m)');
zlabel('Height (m)');
%--------------------------------------------------------------------------