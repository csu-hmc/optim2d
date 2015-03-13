function [R_JRF, L_JRF] = getJointReactionForces(model, x)

% Function to determine the joint reaction forces given the ground reaction
% forces and the forces in the muscles
% Inputs: model:        the model that is used
%         x:            the state, used to calculate the GRF and
%         muscleforces
%         ?? directional vector of each muscle


data = xlsread(model.parameterfile,'','','basic');
LinkMass = data(14:17,2);
g = 9.81;
N = size(x,2);

for i = 1:N
    GRF = gait2d('GRF',x(1:50,i));
    muscleforces = model.gait2d('Muscleforces',x(1:50,i));
    muscleforces = reshape(muscleforces,8,2);
    
    % Get orientation of the femur and tibia
    d = model.gait2d('Stick',x(1:50,i));
    
    OR_femr = atan2(d(2,2)-d(3,2),d(2,1)-d(3,1));
    OR_feml = atan2(d(2,2)-d(7,2),d(2,1)-d(7,1));
    OR_tibr = atan2(d(3,2)-d(4,2),d(3,1)-d(4,1));
    OR_tibl = atan2(d(7,2)-d(8,2),d(7,1)-d(8,1));
    

    % For right leg
    Fresankle = GRF([1,2]) - [0;LinkMass(end)*g];
    Fmusankle = sum(muscleforces(end-2:end,1))*[cos(OR_tibr);sin(OR_tibr)];

    R_JRF(1:2,i) = -Fresankle-Fmusankle;

    Fresknee = -R_JRF(1:2,i) -[0; LinkMass(end-1)*g];
    Fmusknee = sum(muscleforces(3:5,1))*[cos(OR_tibr);sin(OR_tibr)]-sum(muscleforces(7:8,1))*[cos(OR_tibr);sin(OR_tibr)];

    R_JRF(3:4,i) = -Fresknee-Fmusknee;

    Freship = -R_JRF(3:4,i) - [0;LinkMass(end-2)*g];
    Fmuship = sum(muscleforces(1:2,1))*[cos(OR_femr);sin(OR_femr)]-sum(muscleforces([3,5,6],1))*[cos(OR_femr);sin(OR_femr)];

    R_JRF(5:6,i) = -Freship-Fmuship;
    
    % For left leg
    Fresankle = GRF([3,4]) - [0;LinkMass(end)*g];
    Fmusankle = sum(muscleforces(end-2:end,2))*[cos(OR_tibl);sin(OR_tibl)];

    L_JRF(1:2,i) = -Fresankle-Fmusankle;

    Fresknee = -L_JRF(1:2,i) -[0; LinkMass(end-1)*g];
    Fmusknee = sum(muscleforces(3:5,2))*[cos(OR_tibl);sin(OR_tibl)]-sum(muscleforces(7:8,2))*[cos(OR_tibl);sin(OR_tibl)];

    L_JRF(3:4,i) = -Fresknee-Fmusknee;

    Freship = -L_JRF(3:4,i) - [0;LinkMass(end-2)*g];
    Fmuship = sum(muscleforces(1:2,2))*[cos(OR_feml);sin(OR_feml)]-sum(muscleforces([3,5,6],2))*[cos(OR_feml);sin(OR_feml)];

    L_JRF(5:6,i) = -Freship-Fmuship;
end
