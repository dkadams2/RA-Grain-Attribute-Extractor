%Calculate the RA Schmid Factor
function [RASchmidFactor] = RA_SF_calculator(numRAgrains, RAGrainFile)

%Create Stress Direction Matrix (Tensile sample is pulled in the
%x-direction, so it will be [1 0 0]
TensileDir = [1;0;0];
%Create rotation matrix from OpenXY code for each RA grain
for n=1:numRAgrains
    %Convert Euler angles from degrees to radians
    phi1 = pi/180*RAGrainFile(n,2);
    PHI = pi/180*RAGrainFile(n,3);
    phi2 = pi/180*RAGrainFile(n,4);
    
    RAgmat = euler2gmat(phi1,PHI,phi2);
    
    %Create sigma matrix (to put the stress matrix into the crystal frame?)
    sigma = RAgmat*TensileDir;

    %Create FCC Slip Systems Matrix
    FCCSlip = [1, 1, 1, 1, 0, -1;
               1, 1, 1, 1,-1,  0;
               1, 1, 1, 0, 1, -1;
              -1, 1, 1, 1, 1,  0;
              -1, 1, 1, 0, 1, -1;
              -1, 1, 1, 1, 0,  1;
               1,-1, 1, 0, 1,  1;
               1,-1, 1,-1, 0,  1;
               1,-1, 1, 1, 1,  0;
               1, 1,-1, 1, 0,  1;
               1, 1,-1, 0, 1,  1;
               1, 1,-1, 1,-1,  0]; 


    %Solve for angles for Schmid factor
        for p = 1:12
           slipplaneFCC = FCCSlip(p,1:3);
           slipdirFCC = FCCSlip(p,4:6);
           
           thetaFCC = acos(dot(sigma,slipplaneFCC)/(norm(sigma)*norm(slipplaneFCC)));
           lambdaFCC = acos(dot(sigma,slipdirFCC)/(norm(sigma)*norm(slipdirFCC)));
           schmidFCC(p) = cos(thetaFCC)*cos(lambdaFCC);
        end
     RASchmidFactor(n) = max(abs(schmidFCC));
end