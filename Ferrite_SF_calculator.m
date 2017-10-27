%Calculate the RA Schmid Factor
function [FerriteSchmidFactor] = Ferrite_SF_calculator(FerriteGrainFile)

%Create Stress Direction Matrix (Tensile sample is pulled in the
%x-direction, so it will be [1 0 0]
TensileDir = [1;0;0];

%Create rotation matrix from OpenXY code for each RA grain
    %Convert Euler angles from degrees to radians
    phi1 = pi/180*FerriteGrainFile(1,2);
    PHI = pi/180*FerriteGrainFile(1,3);
    phi2 = pi/180*FerriteGrainFile(1,4);
    
    FerriteGmat = euler2gmat(phi1,PHI,phi2);
    
    %Create sigma matrix (to put the stress matrix into the crystal frame?)
    sigma = FerriteGmat*TensileDir;

    %Create BCC Slip Systems Matrix
    BCCSlip = [
    1   1   0  -1   1   1   %1  (110)<-111>
    1   1   0   1  -1   1   %2
    1  -1   0   1   1   1   %3
    1   -1  0   1   1   -1  %4
    1   0   1   1   1   -1  %5
    1   0   1   -1  1   1   %6
    1   0   -1  1   1   1   %7
    1   0   -1  1   -1  1   %8
    0   1   1   1   1   -1  %9
    0   1   1   1   -1  1   %10
    0   1   -1  1   1   1   %11
    0   1   -1  -1  1   1   %12
    1   1   2   1   1   -1  %13
    -1  1   2   1   -1  1   %14
    1   -1  2   -1  1   1   %15
    1   1   -2  1   1   1   %16
    1   2   1   1   -1  1   %17
    -1  2   1   1   1   -1  %18
    1   -2  1   1   1   1   %19
    1   2   -1  -1  1   1   %20
    2   1   1   -1  1   1   %21
    -2  1   1   1   1   1   %22
    2   -1  1   1   1   -1  %23
    2   1   -1  1   -1  1   %24
    1   2   3   1   1   -1  %25
    -1  2   3   1   -1  1   %26
    1   -2  3   -1  1   1   %27
    1   2   -3  1   1   1   %28
    1   3   2   1   -1  1   %29
    -1  3   2   1   1   -1  %30
    1   -3  2   1   1   1   %31
    1   3   -2  -1  1   1   %32
    2   1   3   1   1   -1  %33
    -2  1   3   1   -1  1   %34
    2   -1  3   -1  1   1   %35
    2   1   -3  1   1   1   %36
    2   3   1   1   -1  1   %37
    -2  3   1   1   1   -1  %38
    2   -3  1   1   1   1   %39
    2   3   -1  -1  1   1   %40
    3   1   2   -1  1   1   %41
    -3  1   2   1   1   1   %42
    3   -1  2   1   1   -1  %43
    3   1   -2  1   -1  1   %44
    3   2   1   -1  1   1   %45
    -3  2   1   1   1   1   %46
    3   -2  1   1   1   -1  %47
    3   2   -1  1   -1  1   %48
    ];


    %Solve for angles for Schmid factor
        for p = 1:48
           slipplaneBCC = BCCSlip(p,1:3);
           slipdirBCC = BCCSlip(p,4:6);
           
           thetaBCC = acos(dot(sigma,slipplaneBCC)/(norm(sigma)*norm(slipplaneBCC)));
           lambdaBCC = acos(dot(sigma,slipdirBCC)/(norm(sigma)*norm(slipdirBCC)));
           schmidBCC(p) = cos(thetaBCC)*cos(lambdaBCC);
        end
     FerriteSchmidFactor = max(abs(schmidBCC));
