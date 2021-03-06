%RA Grain Data Compiler
%From: Derrik Adams
%Date: 8/24/17
%This code imports a customized grain file from OpenXY and exports an excel
%file that contains attributes of the RA grains in the scan and their
%surrounding grains.

%Reguired information in grain file: 
%ID, Orientation (in degrees), Position (x,y in microns),Image Quality, Integer identifying phase,
%Area in square microns, Diameter in microns, Aspect Ratio, Major Axis
%Orientation, Neighbors: Count & ID's

%Required Taylor Factor Map: 
%Set Ferrite TF values to the grayscale map
%Set Austenite TF values to the color map

clear all;

%Load Grain File
grainfilename = uigetfile({'*.txt'},'Select Grain File');
grainfile = dlmread(grainfilename);
[grainfilesizeRows, grainfilesizeCols] = size(grainfile);

%Load TF Map File
MapFilename = uigetfile({'*.txt'},'Select Taylor Factor Map File');
TFmapfile = dlmread(MapFilename);
[TFmapsizeRows, TFmapsizeCols] = size(TFmapfile);

%Loop through the grain file to get the coordinates needed and then find
%those coordinates in the map file and assign the TF vector as such
for i=1:grainfilesizeRows
    Xpos = grainfile(i,5); %x-position of the grain 
    Ypos = grainfile(i,6); %y-position of the grain
    for j = 1:TFmapsizeRows
        
        %Need to make ranges on the variation of Xpos and Ypos of 1/2 the
        %step size in um since the map coordinates are by each point (so by
        %the step size) but the grain file coordinates are not on the same
        %step size, they are custom to position but fall in the correct
        %range...
        stepsize = .08; %in um
        range = .5*stepsize;
        if  (TFmapfile(j,1)-range) <= Xpos && Xpos <= (TFmapfile(j,1)+range) && (TFmapfile(j,2)-range) <= Ypos && Ypos <= (TFmapfile(j,2)+range) 
            if TFmapfile(j,3) ~= -1
                TFgrainfile(i) = TFmapfile(j,3);
            else
                TFgrainfile(i) = TFmapfile(j,4);
            end
        end
    end
end

%Append the grainfile to add on the TFgrainfile at end
grainfile(:,grainfilesizeCols+1) = TFgrainfile(1,:);

%Create Variable Names from grainfile
ID = grainfile(:,1);    %Grain ID number
phi1 = grainfile(:,2);  %in deg
PHI = grainfile(:,3);   %in deg
phi2 = grainfile(:,4);  %in deg
Xposition = grainfile(:,5); %in microns
Yposition = grainfile(:,6); %in microns
IQ = grainfile(:,7);    %IQ number
Phase = grainfile(:,8); %1-Austenite; 2-Ferrite
Area = grainfile(:,9);  %in microns
Diameter = grainfile(:,10); %in mircons
GSAR = grainfile(:,11); %Grain Shape Aspect Ratio
MAO = grainfile(:,12);  %Major Axis Orientation in degrees
NumNeighbor = grainfile(:,13); %Number of neighbors

%Extract given data for only austenite grains
index = 1;
for i=1:grainfilesizeRows
    %The phase number will depend on the grain file, Austenite may be 2
    if Phase(i) == 1
        RAgrainLocation(index) = i;
        index = index+1;
    end
end
numRAgrains = length(RAgrainLocation);%defines the number of RA grains identified

%Create a grain file for only the austenite grains
for i = 1:numRAgrains
    RAGrainFileALL(i,:) = grainfile(RAgrainLocation(i),:);
end

for k=1:numRAgrains 
    RAGrainFileCustom(k,1) = grainfile(RAgrainLocation(k),1);  %ID
    RAGrainFileCustom(k,2) = grainfile(RAgrainLocation(k),2);  %phi1
    RAGrainFileCustom(k,3) = grainfile(RAgrainLocation(k),3);  %PHI
    RAGrainFileCustom(k,4) = grainfile(RAgrainLocation(k),4);  %phi2
    RAGrainFileCustom(k,5) = grainfile(RAgrainLocation(k),5);  %X position
    RAGrainFileCustom(k,6) = grainfile(RAgrainLocation(k),6);  %Y position
    RAGrainFileCustom(k,7) = grainfile(RAgrainLocation(k),9);  %Area    
    RAGrainFileCustom(k,8) = grainfile(RAgrainLocation(k),10); %Diameter    
    RAGrainFileCustom(k,9) = grainfile(RAgrainLocation(k),11); %GSAR    
    RAGrainFileCustom(k,10) = grainfile(RAgrainLocation(k),12); %MAO    
    RAGrainFileCustom(k,11) = grainfile(RAgrainLocation(k),13); %NumNeighbor
    RAGrainFileCustom(k,12) = grainfile(RAgrainLocation(k),grainfilesizeCols+1); %TF
end

%Changed from m-2 to m-1 when I added TF value 10/27/17
%14 is the column that neighbor ID's begin at in the original grainfile.
%This adds the neighbor ID values to the end of the RAgrainfilecustom
for m=14:grainfilesizeCols
    RAGrainFileCustom(:,m-1) = RAGrainFileALL(:,m); %Add the neighbor IDs to end
end

%Call Function to compute SCHMID FACTOR for RA grains
RASchmidFactor = RA_SF_calculator(numRAgrains, RAGrainFileCustom);


%% IQ Cutoff GUI Creation
%Have user define IQ value for cutoff between martensite and ferrite.
f = figure('Name','IQ Cutoff Value','Position', [750 500 450 150]);
t = uitable(f, 'ColumnWidth',{150});
x = (max(IQ)-min(IQ))/2; %Default IQ cutoff

IQmin = num2str(min(IQ));
IQmax = num2str(max(IQ));

displaysuggestions = {'Please input the desired IQ cutoff value to distinguish Martensite from Ferrite, where phases identified by OIM with an IQ higher than the input value will be identified as Ferrite and those lower than the input value will be identified as Martensite. There is not a hard cutoff value as IQ values may change from one scan to the next.'};
uilegend = uicontrol('Style','Text');
set(uilegend, 'string', displaysuggestions, 'Position', [25 70 425 75]);

displaymin = {'Minimum IQ for Scan', IQmin};
minlegend = uicontrol('Style','Text');
set(minlegend, 'string', displaymin, 'Position', [210 -40 120 100]);

displaymax = {'Maximum IQ for Scan', IQmax};
minlegend = uicontrol('Style','Text');
set(minlegend, 'string', displaymax, 'Position', [330 -40 120 100]);

cnames = {'Desired IQ Cutoff Value'};
set (t, 'Data', x, 'ColumnEditable', true(1,2), 'ColumnName', cnames);
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

results = menu('Save Changes?', 'OK');
if (results == 1)
    MFcutoff = get(t,'Data');
    close(f)
end
%% RA Grain Neighbor Analysis
%Determine the phases around the RA grain by identifying the surrounding
%grains and their IQ values. Also, calculate their Schmid Factor and Taylor Factors.

%1 Austenite; 2-Ferrite; 3-Martensite

%Create an array of only the column that holds the number of neighbors for
%RA grains
RANumNeighbor = RAGrainFileCustom(:,11);

for n=1:length(RANumNeighbor)
    %Check to see if the RA grain is reported to have 0 neighbors or not.
    %If it has more than 0 neighbors then analyze those neighbors,
    %otherwise set all of the neighbor information to 0 (see else if below)
    if RANumNeighbor(n) ~= 0 
        
        %Variable for loop, will loop through the number of times equal to
        %the number of neighbors the RA grain has.
        for p=1:RANumNeighbor(n)
            %Build a matrix to hold the ID values of the neighboring
            %grains. 12 is the last column number, it holds the TF values,
            %and 13 on hold the neighbor ID values.
            RANeighborID(n,p) = RAGrainFileCustom(n,12+p);
            
            %Loop through the grainfile to find the corresponding ID values
            %and get the info necessary from them.
            for q=1:length(grainfile)
                if RANeighborID(n,p) == ID(q)
                   
                   %Create matrix of initial phase values only according to
                   %the phases scanned for (in our case 1-austenite,
                   %2-ferrite)
                   RANeighborGrainFilePhase(n,p) = Phase(q);
                   
                   %Create a matrix of TF values for neighboring grains
                   RANeighborTF(n,p) = grainfile(q,grainfilesizeCols+1);
                   
                   %Calculate the SF value based on what phase is reported
                   if RANeighborGrainFilePhase(n,p) == 1
                       RANeighborPhaseSchmid(n,p) = RA_SF_calculator(1,grainfile(q,:));
                   else
                       RANeighborPhaseSchmid(n,p) = Ferrite_SF_calculator(grainfile(q,:));
                   end
                   
                   %Store the values of the area of each neighbor grain
                   %reported to be surrounding the RA grain
                   RANeighborGrainArea(n,p) = Area(q);
                   
                   %IF loops to separate the phases based on IQ and Schmid
                   %Factor (Used to use SF to create cutoff between hard
                   %and soft phase, doesn't do that but can be modified to
                   %do so)
                   if Phase(q) == 1 && RANeighborPhaseSchmid(n,p) > .45
                       RANeighborPhase(n,p) = 1; %Austenite
                   end
                   if Phase(q) == 1 && RANeighborPhaseSchmid(n,p) <= .45
                       RANeighborPhase(n,p) = 1; %Austenite
                   end
                   if Phase(q) == 2 && IQ(q)> MFcutoff && RANeighborPhaseSchmid(n,p) > .45
                       RANeighborPhase(n,p) = 2; %Ferrite
                   end
                   if Phase(q) == 2 && IQ(q)> MFcutoff && RANeighborPhaseSchmid(n,p) <= .45
                       RANeighborPhase(n,p) = 2; %Ferrite
                   end
                   if Phase(q) == 2 && IQ(q)<= MFcutoff && RANeighborPhaseSchmid(n,p) > .45
                       RANeighborPhase(n,p) = 3; %Martensite
                   end
                   if Phase(q) == 2 && IQ(q)<= MFcutoff && RANeighborPhaseSchmid(n,p) <= .45
                       RANeighborPhase(n,p) = 3; %Martensite
                   end
                  
                end
            end
        end
    else
        %Set neighbor properties of RA grains reported to have 0 neighbors
        %to all 0's. Helps maintain the correct number of RA grains.
       for p=1:1
            RANeighborID(n,p) = RAGrainFileCustom(n,12+p);
            RANeighborGrainFilePhase(n,p) = 0;   
            RANeighborPhase(n,p) = 0;
            RANeighborPhaseSchmid(n,p) = 0;
       end
    end
end

%% Calculate the average TF of the neighboring grains and subtract it from the RA TF (Knezevik paper)
[RANeighborGrainArea_Row, RANeighborGrainArea_Col] = size(RANeighborGrainArea);

RANeighborAreaSum = zeros(RANeighborGrainArea_Row,1);
RANeighborTFsumtotal = zeros(RANeighborGrainArea_Row,1);
for i=1:size(RAGrainFileCustom)
    for j=1:RANeighborGrainArea_Col
        RANeighborTFSum(i,j) = RANeighborTF(i,j)*RANeighborGrainArea(i,j);
        RANeighborAreaSum(i) = RANeighborGrainArea(i,j)+RANeighborAreaSum(i);
        RANeighborTFsumtotal(i) = RANeighborTFSum(i,j)+RANeighborTFsumtotal(i);
    end
    
    RANeighborTFAverage(i) = RANeighborTFsumtotal(i)/RANeighborAreaSum(i);
    deltaTF(i) = RANeighborTFAverage(i)-RAGrainFileCustom(i,12);
    
end
xaxis = 1:size(RAGrainFileCustom);
scatter(xaxis,deltaTF);
xlabel('RA Grain #');
ylabel('Delta M (Mbar-M)');




%% Neighbor Grain Sorter and Size Determiner/Locator

%Use the sort function to arrange the areas in descending order of size.
RANeighborGrainArea = sort(RANeighborGrainArea,2,'descend');

%For Loop to create vectors for the total size of surrounding neighbors and
%approximate 'common boundary' between the grain and the RA grain
for i=1:RANeighborGrainArea_Row
    %Variable to hold the total area of the 6 largest grains surrounding
    %the RA grain
    RATotalNeighborArea = 0;
    
    %For loop to calculate the total area of only the 6 largest surrounding
    %grains. If you want to look at more, change the bound on j to a larger
    %number.
    for j=1:RANeighborGrainArea_Col
        if j<=RANeighborGrainArea_Col && j<=6
            RATotalNeighborArea = RATotalNeighborArea + RANeighborGrainArea(i,j);
        else 
            break
        end
    end
    
    
    RAGrainTotalNeighborArea(i) = RATotalNeighborArea;
    
    %For loop to calculate the % of each surrounding grain to provide an
    %estimate of how much each grain surrounds the RA grain.
    for k=1:RANeighborGrainArea_Col
        if k<=RANeighborGrainArea_Col && k<=6
            NeighborPercentSurrounding(i,k) = RANeighborGrainArea(i,k)/RAGrainTotalNeighborArea(i);
        else
            break
        end
    end
end

%Neighbor Grain Locator
%Initialize a matrix to hold the values of the sorted ID values
RANeighborIDsSorted = zeros(RANeighborGrainArea_Row,RANeighborGrainArea_Col);

for i=1:RANeighborGrainArea_Row
    for j=1:RANeighborGrainArea_Col
        for k=1:RANeighborGrainArea_Col
            
            %Check to see if the ID is 0, meaning that it has no neighbors
            if RANeighborID(i,k) ~=0
                
                 %Loop through the grainfile to find the grain with the
                 %matching area value
                 if RANeighborGrainArea(i,j) == grainfile(RANeighborID(i,k),9)
                     
                     %Check for duplicate grains being reported
                     for m=1:RANeighborGrainArea_Col
                         if RANeighborID(i,m) ~=0;
                            if RANeighborIDsSorted(i,m) == RANeighborID(i,k)
                                 duplicate = 0;
                                 break
                            else
                                 duplicate = 1;
                            end
                         end
                     end
                     
                     if duplicate == 0
                         break
                     else
                         RANeighborIDsSorted(i,j) = RANeighborID(i,k);
                         RANeighborPhaseSorted(i,j) = RANeighborPhase(i,k);
                         RANeighborTFSorted(i,j) = RANeighborTF(i,k);
                     end      
                 end
            end
        end
    end
end

%Create matrices that hold the x and y positions of the center of the
%neighboring grains
for i=1:RANeighborGrainArea_Row
    for j=1:RANeighborGrainArea_Col
        if j<=RANeighborGrainArea_Col && j<=6
            if RANeighborIDsSorted(i,j) ~= 0
                NeighborPositionX(i,j) = Xposition(RANeighborIDsSorted(i,j));
                NeighborPositionY(i,j) = Yposition(RANeighborIDsSorted(i,j));
            end
        else
            break
        end
    end
end

%Surrounding Phase Percent Calculator
%Identify the phase associated with each number in
%NeighborPercentSurrounding, add all of the common phases numbers, then add
%columns into the output file accordingly
SurroundingRApercent = zeros(RANeighborGrainArea_Row,1);
SurroundingFerritePercent = zeros(RANeighborGrainArea_Row,1);
SurroundingMartPercent = zeros(RANeighborGrainArea_Row,1);
for i=1:RANeighborGrainArea_Row
    for k=1:RANeighborGrainArea_Col
        if k<=RANeighborGrainArea_Col && k<=6
            if RANeighborPhaseSorted(i,k) == 1 
                SurroundingRApercent(i) = NeighborPercentSurrounding(i,k) + SurroundingRApercent(i);
            end
            if RANeighborPhaseSorted(i,k) == 2 
                SurroundingFerritePercent(i) = NeighborPercentSurrounding(i,k) + SurroundingFerritePercent(i);
            end
            if RANeighborPhaseSorted(i,k) == 3
                SurroundingMartPercent(i) = NeighborPercentSurrounding(i,k) + SurroundingMartPercent(i);
            end
        end
    end
end


%Split the area around a grain into 8 Regions and sort the neighboring
%grains into those regions. 
%Regions:
%1-Top Left    2-Top   3-Top Right
%4-Left    GRAIN       5-Right
%6-Bot. Left   7-Bot.  8-Bot. Right
%9-RA grain is surrounded by that grain
for i=1:RANeighborGrainArea_Row
    for j=1:RANeighborGrainArea_Col
        if j<=RANeighborGrainArea_Col && j<=6
            
        %Percent is a fraction of the RA grains diameter. This means that each
        %grain will have different size sectors and it can be changed to
        %provide different location results based on how you'd like to do it. 
        percent = RAGrainFileCustom(i,8)*.25;
        
        %Difference between the x and y positions of the center of the
        %neighbor and the center of the RA grain
        diffX = NeighborPositionX(i,j)-RAGrainFileCustom(i,5);
        diffY = NeighborPositionY(i,j)-RAGrainFileCustom(i,6);
            if NeighborPositionX(i,j) ~=0
                 if diffX < (-1*percent) && diffY < (-1*percent)
                      NeighborLocation(i,j) = 1;
                 end
                 if diffX > (-1*percent) && diffX < percent && diffY < percent
                      NeighborLocation(i,j) = 2;
                 end
                 if diffX > percent && diffY < percent
                      NeighborLocation(i,j) = 3;
                 end
                 if diffX < (-1*percent) && diffY > (-1*percent) && diffY < (percent)
                      NeighborLocation(i,j) = 4;
                 end
                 if diffX > percent && diffY > (-1*percent) && diffY < (percent)
                      NeighborLocation(i,j) = 5;
                 end
                 if diffX < (-1*percent) && diffY > percent
                      NeighborLocation(i,j) = 6;
                 end
                 if diffX > (-1*percent) && diffX < percent && diffY > percent
                      NeighborLocation(i,j) = 7;
                 end
                 if diffX > percent && diffY > percent
                      NeighborLocation(i,j) = 8;
                 end
                 if RAGrainFileCustom(i,11) == 1
                     NeighborLocation(i,j) = 9;
                 end
            end
        else
            break
        end
    end
end
%% Populate the Area Sector vectors
%Using the neighbor location combined with the neighbor phase sorted
%values, populate the vectors to be added into the file RA grain file
Sector1_Phase = zeros(RANeighborGrainArea_Row,1);
Sector2_Phase = zeros(RANeighborGrainArea_Row,1);
Sector3_Phase = zeros(RANeighborGrainArea_Row,1);
Sector4_Phase = zeros(RANeighborGrainArea_Row,1);
Sector5_Phase = zeros(RANeighborGrainArea_Row,1);
Sector6_Phase = zeros(RANeighborGrainArea_Row,1);
Sector7_Phase = zeros(RANeighborGrainArea_Row,1);
Sector8_Phase = zeros(RANeighborGrainArea_Row,1);
Sector1_TF = zeros(RANeighborGrainArea_Row,1);
Sector2_TF = zeros(RANeighborGrainArea_Row,1);
Sector3_TF = zeros(RANeighborGrainArea_Row,1);
Sector4_TF = zeros(RANeighborGrainArea_Row,1);
Sector5_TF = zeros(RANeighborGrainArea_Row,1);
Sector6_TF = zeros(RANeighborGrainArea_Row,1);
Sector7_TF = zeros(RANeighborGrainArea_Row,1);
Sector8_TF = zeros(RANeighborGrainArea_Row,1);

for i=1:RANeighborGrainArea_Row
    for j=1:RANeighborGrainArea_Col
        if j<=RANeighborGrainArea_Col && j<=6
            location = NeighborLocation(i,j);
            switch location
                %The 'if' statements inside the cases are for the times
                %where there are multiple phases assigned to the same
                %vector. Without the 'if' statements, the phase value would
                %correspond to the smallest grain, but it needs to be the
                %largest to maintain correct representation.
                case 1
                   if Sector1_Phase(i) == 0
                    Sector1_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector1_TF(i) = RANeighborTFSorted(i,j);
                   end                     
                case 2
                    if Sector2_Phase(i) == 0
                    Sector2_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector2_TF(i) = RANeighborTFSorted(i,j);
                    end
                case 3
                    if Sector3_Phase(i) == 0
                    Sector3_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector3_TF(i) = RANeighborTFSorted(i,j);
                    end
                case 4
                    if Sector4_Phase(i) == 0
                    Sector4_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector4_TF(i) = RANeighborTFSorted(i,j);
                    end
                case 5
                    if Sector5_Phase(i) == 0
                    Sector5_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector5_TF(i) = RANeighborTFSorted(i,j);
                    end
                case 6
                    if Sector6_Phase(i) == 0
                    Sector6_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector6_TF(i) = RANeighborTFSorted(i,j);
                    end
                case 7
                    if Sector7_Phase(i) == 0
                    Sector7_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector7_TF(i) = RANeighborTFSorted(i,j);
                    end
                case 8
                    if Sector8_Phase(i) == 0
                    Sector8_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector8_TF(i) = RANeighborTFSorted(i,j);
                    end
                case 9
                    %Case 9 says that the grain is surrounded by its
                    %neighbor, so assign all sectors to the same phase and
                    %TF value.
                    Sector1_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector2_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector3_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector4_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector5_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector6_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector7_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector8_Phase(i) = RANeighborPhaseSorted(i,j);
                    Sector1_TF(i) = RANeighborTFSorted(i,j);
                    Sector2_TF(i) = RANeighborTFSorted(i,j);
                    Sector3_TF(i) = RANeighborTFSorted(i,j);
                    Sector4_TF(i) = RANeighborTFSorted(i,j);
                    Sector5_TF(i) = RANeighborTFSorted(i,j);
                    Sector6_TF(i) = RANeighborTFSorted(i,j);
                    Sector7_TF(i) = RANeighborTFSorted(i,j);
                    Sector8_TF(i) = RANeighborTFSorted(i,j);
            end           
         else
            break
        end
    end
end
%% Blank Sector Filler
%Fill in the blank sector spaces with the phase of the largest neighboring
%sector
for i=1:RANeighborGrainArea_Row
    %Make them all 0 for each RA grain so it doesn't carry on to the next
    %loop
    size1 = 0;
    size2 = 0;
    size3 = 0;
    size4 = 0;
    size5 = 0;
    size6 = 0;
    size7 = 0;
    size8 = 0;
    
    for j=1:RANeighborGrainArea_Col
        if j<=RANeighborGrainArea_Col && j<=6
            location = NeighborLocation(i,j);
            
            switch location
                %The 'if' statements inside the cases are for the cases
                %where there are multiple phases assigned to the same
                %vector. Without the 'if' statements, the size value would
                %be assigned to the smallest grain, but it needs to be the
                %largest to maintain correct representation.
                case 1
                    if size1 == 0
                    size1 = RANeighborGrainArea(i,j);
                    end
                case 2
                    if size2 == 0
                    size2 = RANeighborGrainArea(i,j);
                    end
                case 3
                    if size3 == 0
                    size3 = RANeighborGrainArea(i,j);
                    end
                case 4
                    if size4 == 0
                    size4 = RANeighborGrainArea(i,j);
                    end
                case 5
                    if size5 == 0
                    size5 = RANeighborGrainArea(i,j);
                    end
                case 6
                    if size6 == 0
                    size6 = RANeighborGrainArea(i,j);
                    end
                case 7
                    if size7 == 0
                    size7 = RANeighborGrainArea(i,j);
                    end
                case 8
                    if size8 == 0
                    size8 = RANeighborGrainArea(i,j);
                    end
            end
        end 
    end
    %This is rather tedious code, but basically it will start at sector 1
    %and work its way through each sector to see if they have a phase
    %assigned to them. If not, it will compare the two immediate
    %neighboring sectors and assign the phase of the larger of the two. If
    %those are both empty it will move on to the next two sectors and so on
    %until it finds a phase. 
     if Sector1_Phase(i) == 0
         if size2 ~=0 || size4 ~=0
             if size2 > size4 
                 Sector1_Phase(i) = Sector2_Phase(i);
                 Sector1_TF(i) = Sector2_TF(i);
             else
                 Sector1_Phase(i) = Sector4_Phase(i);
                 Sector1_TF(i) = Sector4_TF(i);
             end
         end
             if Sector1_Phase(i) == 0 
                 if size3 ~=0 || size6 ~=0
                    if size6 > size3 
                        Sector1_Phase(i) = Sector6_Phase(i);
                        Sector1_TF(i) = Sector6_TF(i);
                    else
                        Sector1_Phase(i) = Sector3_Phase(i);
                        Sector1_TF(i) = Sector3_TF(i);
                    end
                 end
              if Sector1_Phase(i) == 0 
                        if size3 ~=7 || size5 ~=0
                            if size7 > size5
                                Sector1_Phase(i) = Sector7_Phase(i);
                                Sector1_TF(i) = Sector7_TF(i);
                            else
                                Sector1_Phase(i) = Sector5_Phase(i);
                                Sector1_TF(i) = Sector5_TF(i);
                            end
                        else
                            Sector1_Phase(i) = Sector8_Phase(i);
                            Sector1_TF(i) = Sector8_TF(i);
                        end
               end
             end
     end
      if Sector2_Phase(i) == 0
          if size3 ~=0 || size1 ~=0
             if size3 > size1  
                 Sector2_Phase(i) = Sector3_Phase(i);
                 Sector2_TF(i) = Sector3_TF(i);
             else
                 Sector2_Phase(i) = Sector1_Phase(i);
                 Sector2_TF(i) = Sector1_TF(i);
             end
          end
             if Sector2_Phase(i) == 0
                 if size4 ~=0 || size5 ~=0
                        if size4 > size5 
                            Sector2_Phase(i) = Sector4_Phase(i);
                            Sector2_TF(i) = Sector4_TF(i);
                        else
                            Sector2_Phase(i) = Sector5_Phase(i);
                            Sector2_TF(i) = Sector5_TF(i);
                        end
                 end
                 if Sector2_Phase(i) == 0
                        if size6 ~=0 || size8 ~=0
                         
                            if size6 > size8
                                Sector2_Phase(i) = Sector6_Phase(i);
                                Sector2_TF(i) = Sector6_TF(i);
                            else
                                Sector2_Phase(i) = Sector8_Phase(i);
                                Sector2_TF(i) = Sector8_TF(i);
                            end
                        else
                          Sector2_Phase(i) = Sector7_Phase(i);
                          Sector2_TF(i) = Sector7_TF(i);
                        end
                 end
             end
      end
       if Sector3_Phase(i) == 0
           if size2 ~=0 || size5 ~=0
             if size2 > size5
                 Sector3_Phase(i) = Sector2_Phase(i);
                 Sector3_TF(i) = Sector2_TF(i);
             else
                 Sector3_Phase(i) = Sector5_Phase(i);
                 Sector3_TF(i) = Sector5_TF(i);
             end
           end
           if Sector3_Phase(i) == 0
                if size1 ~=0 || size8 ~=0
                    if size1 > size8
                        Sector3_Phase(i) = Sector1_Phase(i);
                        Sector3_TF(i) = Sector1_TF(i);
                    else
                        Sector3_Phase(i) = Sector8_Phase(i);
                        Sector3_TF(i) = Sector8_TF(i);
                    end
                end
           end
                if Sector3_Phase(i) == 0
                    if size4 ~=0 || size7 ~=0
                        if size4 > size7 
                            Sector3_Phase(i) = Sector4_Phase(i);
                            Sector3_TF(i) = Sector4_TF(i);
                        else
                            Sector3_Phase(i) = Sector7_Phase(i);
                            Sector3_TF(i) = Sector7_TF(i);
                        end
                    else
                       Sector3_Phase(i) = Sector6_Phase(i);
                       Sector3_TF(i) = Sector6_TF(i);
                    end
                end
      end
       if Sector4_Phase(i) == 0
           if size6 ~=0 || size1 ~=0
             if size6 > size1
                 Sector4_Phase(i) = Sector6_Phase(i);
                 Sector4_TF(i) = Sector6_TF(i);
             else
                 Sector4_Phase(i) = Sector1_Phase(i);
                 Sector4_TF(i) = Sector1_TF(i);
             end
           end
             if Sector4_Phase(i) == 0
                 if size2 ~=0 || size7 ~=0
                      if size2 > size7 
                          Sector4_Phase(i) = Sector2_Phase(i);
                          Sector4_TF(i) = Sector2_TF(i);
                      else
                          Sector4_Phase(i) = Sector7_Phase(i);
                          Sector4_TF(i) = Sector7_TF(i);
                      end
                 end
             end
                if Sector4_Phase(i) == 0
                    if size3 ~=0 || size8 ~=0
                        if size3 > size8
                            Sector4_Phase(i) = Sector3_Phase(i);
                            Sector4_TF(i) = Sector3_TF(i);
                        else
                            Sector4_Phase(i) = Sector8_Phase(i);
                            Sector4_TF(i) = Sector8_TF(i);
                        end
                    else
                       Sector4_Phase(i) = Sector5_Phase(i);
                       Sector4_TF(i) = Sector5_TF(i);
                    end
                end
            
        end
      if Sector5_Phase(i) == 0
          if size3 ~=0 || size8 ~=0
             if size3 > size8
                 Sector5_Phase(i) = Sector3_Phase(i);
                 Sector5_TF(i) = Sector3_TF(i);
             else
                 Sector5_Phase(i) = Sector8_Phase(i);
                 Sector5_TF(i) = Sector8_TF(i);
             end
          end
             if Sector5_Phase(i) == 0
                 if size2 ~=0 || size7 ~=0
                      if size2 > size7
                         Sector5_Phase(i) = Sector2_Phase(i);
                         Sector5_TF(i) = Sector2_TF(i);
                      else
                         Sector5_Phase(i) = Sector7_Phase(i);
                         Sector5_TF(i) = Sector7_TF(i);
                      end
                 end
             end
                if Sector5_Phase(i) == 0
                    if size1 ~=0 || size6 ~=0
                         if size1 > size6
                             Sector5_Phase(i) = Sector1_Phase(i);
                             Sector5_TF(i) = Sector1_TF(i);
                         else
                             Sector5_Phase(i) = Sector6_Phase(i);
                             Sector5_TF(i) = Sector6_TF(i);
                         end
                     else
                       Sector5_Phase(i) = Sector4_Phase(i);
                       Sector5_TF(i) = Sector4_TF(i);
                    end
                end
      end
       if Sector6_Phase(i) == 0
           if size4 ~=0 || size7 ~=0
             if size4 > size7
                 Sector6_Phase(i) = Sector4_Phase(i);
                 Sector6_TF(i) = Sector4_TF(i);
             else
                 Sector6_Phase(i) = Sector7_Phase(i);
                 Sector6_TF(i) = Sector7_TF(i);
             end
           end
             if Sector6_Phase(i) == 0
                 if size1 ~=0 || size8 ~=0
                        if size1 > size8
                            Sector6_Phase(i) = Sector1_Phase(i);
                            Sector6_TF(i) = Sector1_TF(i);
                        else
                            Sector6_Phase(i) = Sector8_Phase(i);
                            Sector6_TF(i) = Sector8_TF(i);
                        end
                 end
             end
                if Sector6_Phase(i) == 0
                    if size2 ~=0 || size5 ~=0
                        if size2 > size5
                            Sector6_Phase(i) = Sector2_Phase(i);
                            Sector6_TF(i) = Sector2_TF(i);
                        else
                            Sector6_Phase(i) = Sector5_Phase(i);
                            Sector6_TF(i) = Sector5_TF(i);
                        end
                    else
                       Sector6_Phase(i) = Sector3_Phase(i);
                       Sector6_TF(i) = Sector3_TF(i);
                    end
                end
            
       end
       if Sector7_Phase(i) == 0
          if size6 ~=0 || size8 ~=0
             if size6 > size8
                 Sector7_Phase(i) = Sector6_Phase(i);
                 Sector7_TF(i) = Sector6_TF(i);
             else
                 Sector7_Phase(i) = Sector8_Phase(i);
                 Sector7_TF(i) = Sector8_TF(i);
             end
          end
             if Sector7_Phase(i) == 0
                 if size4 ~=0 || size5 ~=0
                     if size4 > size5
                         Sector7_Phase(i) = Sector4_Phase(i);
                         Sector7_TF(i) = Sector4_TF(i);
                     else
                         Sector7_Phase(i) = Sector5_Phase(i);
                         Sector7_TF(i) = Sector5_TF(i);
                     end
                 end
             end
                if Sector7_Phase(i) == 0
                    if size1 ~= 0 || size3 ~=0
                        if size1 > size3
                             Sector7_Phase(i) = Sector1_Phase(i);
                             Sector7_TF(i) = Sector1_TF(i);
                        else
                            Sector7_Phase(i) = Sector3_Phase(i);
                            Sector7_TF(i) = Sector3_TF(i);
                        end
                    else
                        Sector7_Phase(i) = Sector2_Phase(i);
                        Sector7_TF(i) = Sector2_TF(i);
                    end
                end
            
      end
      if Sector8_Phase(i) == 0
          if size5 ~=0 || size7 ~=0
             if size5 > size7
                 Sector8_Phase(i) = Sector5_Phase(i);
                 Sector8_TF(i) = Sector5_TF(i);
             else
                 Sector8_Phase(i) = Sector7_Phase(i);
                 Sector8_TF(i) = Sector7_TF(i);
             end
          end
             if Sector8_Phase(i) == 0
                 if size3 ~=0 || size6 ~=0
                      if size3 > size6
                         Sector8_Phase(i) = Sector3_Phase(i);
                         Sector8_TF(i) = Sector3_TF(i);
                      else
                         Sector8_Phase(i) = Sector6_Phase(i);
                         Sector8_TF(i) = Sector6_TF(i);
                      end
                 end
             end
                if Sector8_Phase(i) == 0
                   if size2 ~=0 || size4 ~=0
                        if size2 > size4
                           Sector8_Phase(i) = Sector2_Phase(i);
                           Sector8_TF(i) = Sector2_TF(i);
                        else
                           Sector8_Phase(i) = Sector4_Phase(i);
                           Sector8_TF(i) = Sector4_TF(i);
                        end
                   else
                       Sector8_Phase(i) = Sector1_Phase(i);
                       Sector8_TF(i) = Sector1_TF(i);
                   end
                end
      end
              
end

%% Add Strain Data from DIC
%load the strain data for the strain being analyzed
strain_grain_file1 = uigetfile({'*.txt'}, 'Select Strain File of Grain Type 1');
grain_file1_strain = dlmread(strain_grain_file1);
point_strains = grain_file1_strain(:,11);
point_ID = grain_file1_strain(:,9);
point_X = grain_file1_strain(:,4);
point_Y = grain_file1_strain(:,5);

norm_point_strains = grain_file1_strain(:,12);

%Get the average normalized strain for an RA grain and its neihbors to see
%which RA grains have higher localized strain. This theoretically should
%alleviate the issue of having to use absolute strain values which vary
%from one scan to the next

%Get the average normalized strain in the RA grains
for i = 1:numRAgrains
    norm_strain_total = 0;
    num_strain_points = 0;
    for j = 1:size(norm_point_strains)
        if point_ID(j) == RAGrainFileALL(i)
            norm_strain_total = norm_strain_total+norm_point_strains(j);
            num_strain_points = num_strain_points + 1;
        end
    end
    RAGrainNormalizedStrainAve(i) = norm_strain_total/num_strain_points;
end

%Get normalized average strain of neighbors around the RA grain
Neighbor_norm_strain_total = 0;
strain_count = 0;
for i = 1:numRAgrains
       Neighbor_norm_strain_total = 0;
       strain_count = 0;
    if RANumNeighbor(i) ~= 0 
    for j = 1:RANumNeighbor(i)
        for k=1:size(norm_point_strains)
                if RANeighborID(i,j) == grain_file1_strain(k,9)
                   Neighbor_norm_strain_total = grain_file1_strain(k,12) + Neighbor_norm_strain_total;
                   strain_count = strain_count+1;
                end
        end
    end
        RANeighborNormalizedStrainAve(i) = Neighbor_norm_strain_total/strain_count;
    else
        RANeighborNormalizedStrainAve(i) = 0;
    end
end

%Just get average absolute strain in RA grain for now
for i = 1:numRAgrains
    strain_total = 0;
    num_strain_points = 0;
    for j = 1:size(point_strains)
        if point_ID(j) == RAGrainFileALL(i)
            strain_total = strain_total+point_strains(j);
            num_strain_points = num_strain_points +1;
        end
    end
    RAGrainStrainAve(i) = strain_total/num_strain_points;
end

%Get average absolute strain of neighbors around the RA grain
Neighbor_strain_total = 0;
strain_count = 0;
for i = 1:numRAgrains
      Neighbor_strain_total = 0;
       strain_count = 0;
    if RANumNeighbor(i) ~= 0 
    for j = 1:RANumNeighbor(i)
        for k=1:size(point_strains)
                if RANeighborID(i,j) == grain_file1_strain(k,9)
                   Neighbor_strain_total = grain_file1_strain(k,11) + Neighbor_strain_total;
                   strain_count = strain_count+1;
                end
        end
    end
        RANeighborStrainAve(i) = Neighbor_strain_total/strain_count;     
    else
        RANeighborStrainAve(i) = 0;
    end
end

figure(8)
scatter(1:i,RAGrainStrainAve)
hold on
scatter(1:i,RANeighborStrainAve)
title('RA DIC Strain Averages')
xlabel('RA Grain')
ylabel('Strain')
legend('RA Grain Average Strain','Surrounding Neighbor Strain Average')

figure(9)
scatter(1:i,RAGrainNormalizedStrainAve)
hold on
scatter(1:i,RANeighborNormalizedStrainAve)
title('RA DIC Normalized Strain Averages')
xlabel('RA Grain')
ylabel('Strain')
legend('RA Grain Normalized Average Strain','Surrounding Neighbor Normalized Strain Average')
%% Calculate the average GND Content in the RA grains
% Load the GND file 
GND_file = uigetfile({'*.txt'}, 'Select GND File');
GND_mat = dlmread(GND_file);
GND_x_point = GND_mat(:,1);
GND_y_point = GND_mat(:,2);
GND_value = GND_mat(:,4);

%Add grain ID onto the GND file
for i=1:size(GND_mat)
    GND_mat(i,5) = point_ID(i);
end

%Get the average GND content in the RA grains
for i = 1:numRAgrains
    GND_total = 0;
    num_GND_points = 0;
    for j = 1:size(GND_mat)
        if GND_mat(j,5) == RAGrainFileALL(i)
            GND_total = GND_total+GND_value(j);
            num_GND_points = num_GND_points + 1;
        end
    end
    RAGrainGNDAve(i) = GND_total/num_GND_points;
end

        
%% Final RA File Creation

for q=1:numRAgrains
    FinalRAFile(q,1) = grainfile(RAgrainLocation(q),1);  %ID
    FinalRAFile(q,2) = grainfile(RAgrainLocation(q),9);  %Area
    FinalRAFile(q,3) = grainfile(RAgrainLocation(q),10); %Diameter
    FinalRAFile(q,4) = grainfile(RAgrainLocation(q),11); %GSAR
    FinalRAFile(q,5) = grainfile(RAgrainLocation(q),12); %MAO
    FinalRAFile(q,6) = RASchmidFactor(q);                %RA Schmid Factor
    FinalRAFile(q,7) = grainfile(RAgrainLocation(q),grainfilesizeCols+1); %RA Taylor Factor
    FinalRAFile(q,8) = grainfile(RAgrainLocation(q),13); %NumNeighbor
    FinalRAFile(q,9) = Sector1_Phase(q);
    FinalRAFile(q,10) = Sector2_Phase(q);
    FinalRAFile(q,11) = Sector3_Phase(q);
    FinalRAFile(q,12) = Sector4_Phase(q);
    FinalRAFile(q,13) = Sector5_Phase(q);
    FinalRAFile(q,14) = Sector6_Phase(q);
    FinalRAFile(q,15) = Sector7_Phase(q);
    FinalRAFile(q,16) = Sector8_Phase(q);
    FinalRAFile(q,17) = SurroundingFerritePercent(q);
    FinalRAFile(q,18) = SurroundingMartPercent(q);
    FinalRAFile(q,19) = SurroundingRApercent(q);
    FinalRAFile(q,20) = Sector1_TF(q);
    FinalRAFile(q,21) = Sector2_TF(q);
    FinalRAFile(q,22) = Sector3_TF(q);
    FinalRAFile(q,23) = Sector4_TF(q);
    FinalRAFile(q,24) = Sector5_TF(q);
    FinalRAFile(q,25) = Sector6_TF(q);
    FinalRAFile(q,26) = Sector7_TF(q);
    FinalRAFile(q,27) = Sector8_TF(q);
    FinalRAFile(q,28) = RANeighborTFAverage(q);
    FinalRAFile(q,29) = deltaTF(q);
    
end

Headers = {'GrainID', 'Area', 'Diameter', 'GSAR', 'MAO', 'Schmid_Factor','Taylor_Factor',...
           'Num_Neighbors','Sector1_Phase','Sector2_Phase','Sector3_Phase','Sector4_Phase','Sector5_Phase','Sector6_Phase',...
           'Sector7_Phase','Sector8_Phase', 'Surrounding_Ferrite', 'Surrounding_Martensite','Surrounding_RA',...
           'Sector1_TF','Sector2_TF','Sector3_TF','Sector4_TF','Sector5_TF','Sector6_TF','Sector7_TF','Sector8_TF',...
           'Neighbor_Weighted_Average_TF','TF_Neighbors_Grain_Diff'};

%Create a table with the FinalRAFile info that has headers to describe
%column data.
T = array2table(FinalRAFile,'VariableNames',Headers);

%Delete all rows (or RA grains) that have 0 neighbors
toDelete = T.Num_Neighbors<1; 
T(toDelete,:) = [];

%Create CSV file for RA data by prompting user for filename. File will be
%saved in the folder where this script is located regardless of the path created when
%making the filename.
savequery = menu('Save to a CSV File?', 'Yes','No');
if (savequery == 1)
   [filename,path] = uiputfile('*.csv','Save File As');
    writetable(T, filename);
end



    

