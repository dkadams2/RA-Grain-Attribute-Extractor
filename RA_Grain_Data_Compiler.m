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

clear all;

%Load File
filename = uigetfile({'*.txt'},'Select Grain File');
grainfile = dlmread(filename);
[sizeRows, sizeCols] = size(grainfile);

%Create Variable Names
ID = grainfile(:,1);
phi1 = grainfile(:,2);  %in deg
PHI = grainfile(:,3);   %in deg
phi2 = grainfile(:,4);  %in deg
Xposition = grainfile(:,5); %in microns
Yposition = grainfile(:,6); %in microns
IQ = grainfile(:,7);
Phase = grainfile(:,8); %1-Austenite; 2-Ferrite
Area = grainfile(:,9);  %in microns
Diameter = grainfile(:,10); %in mircons
GSAR = grainfile(:,11);
MAO = grainfile(:,12);  %in degrees
NumNeighbor = grainfile(:,13);

%Extract given data for austenite grains
index = 1;
for i=1:sizeRows
    if Phase(i) == 1
        RAgrainLocation(index) = i;
        index = index+1;
    end
end
numRAgrains = length(RAgrainLocation);%defines the number of RA grains identified

%Create a grain file for only the austenite grains
for j = 1:numRAgrains
    RAGrainFileALL(j,:) = grainfile(RAgrainLocation(j),:);
end

for k=1:numRAgrains 
    RAGrainFile(k,1) = grainfile(RAgrainLocation(k),1);  %ID
    RAGrainFile(k,2) = grainfile(RAgrainLocation(k),2);  %phi1
    RAGrainFile(k,3) = grainfile(RAgrainLocation(k),3);  %PHI
    RAGrainFile(k,4) = grainfile(RAgrainLocation(k),4);  %phi2
    RAGrainFile(k,5) = grainfile(RAgrainLocation(k),5);  %X position
    RAGrainFile(k,6) = grainfile(RAgrainLocation(k),6);  %Y position
    RAGrainFile(k,7) = grainfile(RAgrainLocation(k),9);  %Area    
    RAGrainFile(k,8) = grainfile(RAgrainLocation(k),10); %Diameter    
    RAGrainFile(k,9) = grainfile(RAgrainLocation(k),11); %GSAR    
    RAGrainFile(k,10) = grainfile(RAgrainLocation(k),12); %MAO    
    RAGrainFile(k,11) = grainfile(RAgrainLocation(k),13); %NumNeighbor    
end

for m=14:sizeCols
    RAGrainFile(:,m-2) = RAGrainFileALL(:,m); %Add the neighbor IDs to end
end

[numRARows, numRACols] = size(RAGrainFile);

%Call Function to compute SCHMID FACTOR for RA
RASchmidFactor = RA_SF_calculator(numRAgrains, RAGrainFile);



%% IQ Cutoff GUI Creation
%Have user define IQ value for cutoff between martensite and ferrite.
f = figure('Name','IQ Cutoff Value','Position', [750 500 450 150]);
t = uitable(f, 'ColumnWidth',{150});
x = (max(IQ)-min(IQ))/2; %Default IQ cutoff

IQmin = num2str(min(IQ));
IQmax = num2str(max(IQ));

displaysuggestions = {'Please input the desired IQ cutoff value to distinguish Martensite from Ferrite, where phases identified by OIM with an IQ higher than the input value will be identified as Ferrite and those lower than the input value will be identified as Martensite. There is not a hard cutoff value as IQ values may change from one scan to the next, but for more recent scans, an IQ cutoff value around 2e6 - 2.5e6 seems reasonable.'};
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
%grains and their IQ values. Also, calculate their Schmid Factor.

%1-Soft Austenite; 2-Hard Austenite; 3-Soft Ferrite; 4-Hard Ferrite; 5-Soft Martensite; 6-Hard Martensite

RANumNeighbor = RAGrainFile(:,11);

for n=1:length(RANumNeighbor)
    if RANumNeighbor(n) ~= 0 
        for p=1:RANumNeighbor(n)
            RANeighbor(n,p) = RAGrainFile(n,11+p);
            for q=1:length(grainfile)
                if RANeighbor(n,p) == ID(q)
                   RANeighborInitialPhase(n,p) = Phase(q);
                   
                   
                   if RANeighborInitialPhase(n,p) == 1
                       RANeighborPhaseSchmid(n,p) = RA_SF_calculator(1,grainfile(q,:));
                   else
                       RANeighborPhaseSchmid(n,p) = Ferrite_SF_calculator(grainfile(q,:));
                   end
                   NeighborArea(n,p) = Area(q);
                   
                   %IF loops to separate the phases based on IQ and Schmid
                   %Factor
                   if Phase(q) == 1 && RANeighborPhaseSchmid(n,p) > .45
                       RANeighborPhase(n,p) = 1; %Soft Austenite
                   end
                   if Phase(q) == 1 && RANeighborPhaseSchmid(n,p) <= .45
                       RANeighborPhase(n,p) = 2; %Hard Austenite
                   end
                   if Phase(q) == 2 && IQ(q)> MFcutoff && RANeighborPhaseSchmid(n,p) > .45
                       RANeighborPhase(n,p) = 3; %Soft Ferrite
                   end
                   if Phase(q) == 2 && IQ(q)> MFcutoff && RANeighborPhaseSchmid(n,p) <= .45
                       RANeighborPhase(n,p) = 4; %Hard Austenite
                   end
                   if Phase(q) == 2 && IQ(q)<= MFcutoff && RANeighborPhaseSchmid(n,p) > .45
                       RANeighborPhase(n,p) = 5; %Soft Martensite
                   end
                   if Phase(q) == 2 && IQ(q)<= MFcutoff && RANeighborPhaseSchmid(n,p) <= .45
                       RANeighborPhase(n,p) = 6; %Hard Martensite
                   end
                  
                end
            end
        end
    else
       for p=1:1
            RANeighbor(n,p) = RAGrainFile(n,11+p);
            RANeighborInitialPhase(n,p) = 0;   
            RANeighborPhase(n,p) = 0;
            RANeighborPhaseSchmid(n,p) = 0;
       end
    end
end
%% Neighbor Grain Sorter and Size Determiner/Locator
RAx = RAGrainFile(:,5);
RAy = RAGrainFile(:,6);
[row, col] = size(NeighborArea);

%Use the sort function to arrange the areas in descending order of size.
NeighborArea = sort(NeighborArea,2,'descend');

%For Loop to create vectors for the total size of surrounding neighbors and
%approximate 'common boundary' between the grain and the RA grain
for i=1:row
    NeighborTotalArea = 0;
    %For loop to calculate the total area of only the 6 largest surrounding
    %grains. If you want to look at more, change the bound on j to a larger
    %number.
    for j=1:col
        if j<=col && j<=6
            NeighborTotalArea = NeighborTotalArea + NeighborArea(i,j);
        else 
            break
        end
    end
    GrainNeighborArea(i) = NeighborTotalArea;
    
    %For loop to calculate the % of each surrounding grain to provide an
    %estimate of how much each grain surrounds the RA grain.
    for k=1:col
        if k<=col && k<=6
            NeighborPercentSurrounding(i,k) = NeighborArea(i,k)/GrainNeighborArea(i);
        else
            break
        end
    end
    %NeighborPercentSurrounding = NeighborArea(i,j)/NeighborTotalArea;
end



% Neighbor Grain Locator

RANeighborIndiceSorted = zeros(row,col);

for i=1:row
    for j=1:col
        for k=1:col
            if RANeighbor(i,k) ~=0
                 if NeighborArea(i,j) == grainfile(RANeighbor(i,k),9)
                     
                     for m=1:col
                         if RANeighbor(i,m) ~=0;
                            if RANeighborIndiceSorted(i,m) == RANeighbor(i,k)
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
                         RANeighborIndiceSorted(i,j) = RANeighbor(i,k);
                         RANeighborPhaseSorted(i,j) = RANeighborPhase(i,k);
                     end      
                 end
            end
        end
    end
end


for i=1:row
    for j=1:col
        if j<=col && j<=6
            if RANeighborIndiceSorted(i,j) ~= 0
                NeighborPositionX(i,j) = Xposition(RANeighborIndiceSorted(i,j));
                NeighborPositionY(i,j) = Yposition(RANeighborIndiceSorted(i,j));
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
SurroundingRApercent = zeros(row,1);
SurroundingFerritePercent = zeros(row,1);
SurroundingMartPercent = zeros(row,1);
for i=1:row
    for k=1:col
        if k<=col && k<=6
            if RANeighborPhaseSorted(i,k) == 1 || RANeighborPhaseSorted(i,k) == 2
                SurroundingRApercent(i) = NeighborPercentSurrounding(i,k) + SurroundingRApercent(i);
            end
            if RANeighborPhaseSorted(i,k) == 3 || RANeighborPhaseSorted(i,k) == 4
                SurroundingFerritePercent(i) = NeighborPercentSurrounding(i,k) + SurroundingFerritePercent(i);
            end
            if RANeighborPhaseSorted(i,k) == 5 || RANeighborPhaseSorted(i,k) == 6
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
for i=1:row
    for j=1:col
        if j<=col && j<=6
        %The percent value can be changed to provide different location
        %results based on how you'd like to do it. 
        
        percent = RAGrainFile(i,8)*.25;
        %percent = .1;
        diffX = NeighborPositionX(i,j)-RAGrainFile(i,5);
        diffY = NeighborPositionY(i,j)-RAGrainFile(i,6);
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
                 if RAGrainFile(i,11) == 1
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
Sector1 = zeros(row,1);
Sector2 = zeros(row,1);
Sector3 = zeros(row,1);
Sector4 = zeros(row,1);
Sector5 = zeros(row,1);
Sector6 = zeros(row,1);
Sector7 = zeros(row,1);
Sector8 = zeros(row,1);
for i=1:row
    for j=1:col
        if j<=col && j<=6
            location = NeighborLocation(i,j);
            switch location
                %The 'if' statements inside the cases are for the times
                %where there are multiple phases assigned to the same
                %vector. Without the 'if' statements, the phase value would
                %correspond to the smallest grain, but it needs to be the
                %largest to maintain correct representation.
                case 1
                   if Sector1(i) == 0
                    Sector1(i) = RANeighborPhaseSorted(i,j);
                   end                     
                case 2
                    if Sector2(i) == 0
                    Sector2(i) = RANeighborPhaseSorted(i,j);
                    end
                case 3
                    if Sector3(i) == 0
                    Sector3(i) = RANeighborPhaseSorted(i,j);
                    end
                case 4
                    if Sector4(i) == 0
                    Sector4(i) = RANeighborPhaseSorted(i,j);
                    end
                case 5
                    if Sector5(i) == 0
                    Sector5(i) = RANeighborPhaseSorted(i,j);
                    end
                case 6
                    if Sector6(i) == 0
                    Sector6(i) = RANeighborPhaseSorted(i,j);
                    end
                case 7
                    if Sector7(i) == 0
                    Sector7(i) = RANeighborPhaseSorted(i,j);
                    end
                case 8
                    if Sector8(i) == 0
                    Sector8(i) = RANeighborPhaseSorted(i,j);
                    end
                case 9
                    Sector1(i) = RANeighborPhaseSorted(i,j);
                    Sector2(i) = RANeighborPhaseSorted(i,j);
                    Sector3(i) = RANeighborPhaseSorted(i,j);
                    Sector4(i) = RANeighborPhaseSorted(i,j);
                    Sector5(i) = RANeighborPhaseSorted(i,j);
                    Sector6(i) = RANeighborPhaseSorted(i,j);
                    Sector7(i) = RANeighborPhaseSorted(i,j);
                    Sector8(i) = RANeighborPhaseSorted(i,j);
            end           
         else
            break
        end
    end
end
%% Blank Sector Filler
%Fill in the blank sector spaces with the phase of the largest neighboring
%sector
for i=1:row
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
    
    for j=1:col
        if j<=col && j<=6
            location = NeighborLocation(i,j);
            
            switch location
                %The 'if' statements inside the cases are for the cases
                %where there are multiple phases assigned to the same
                %vector. Without the 'if' statements, the size value would
                %be assigned to the smallest grain, but it needs to be the
                %largest to maintain correct representation.
                case 1
                    if size1 == 0
                    size1 = NeighborArea(i,j);
                    end
                case 2
                    if size2 == 0
                    size2 = NeighborArea(i,j);
                    end
                case 3
                    if size3 == 0
                    size3 = NeighborArea(i,j);
                    end
                case 4
                    if size4 == 0
                    size4 = NeighborArea(i,j);
                    end
                case 5
                    if size5 == 0
                    size5 = NeighborArea(i,j);
                    end
                case 6
                    if size6 == 0
                    size6 = NeighborArea(i,j);
                    end
                case 7
                    if size7 == 0
                    size7 = NeighborArea(i,j);
                    end
                case 8
                    if size8 == 0
                    size8 = NeighborArea(i,j);
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
     if Sector1(i) == 0
         if size2 ~=0 || size4 ~=0
             if size2 > size4 
                 Sector1(i) = Sector2(i);
             else
                 Sector1(i) = Sector4(i);
             end
         end
             if Sector1(i) == 0 
                 if size3 ~=0 || size6 ~=0
                    if size6 > size3 
                        Sector1(i) = Sector6(i);
                    else
                        Sector1(i) = Sector3(i);
                    end
                 end
              if Sector1(i) == 0 
                        if size3 ~=7 || size5 ~=0
                            if size7 > size5
                                Sector1(i) = Sector7(i);
                            else
                                Sector1(i) = Sector5(i);
                            end
                        else
                            Sector1(i) = Sector8(i);
                        end
               end
             end
     end
      if Sector2(i) == 0
          if size3 ~=0 || size1 ~=0
             if size3 > size1  
                 Sector2(i) = Sector3(i);
             else
                 Sector2(i) = Sector1(i);
             end
          end
             if Sector2(i) == 0
                 if size4 ~=0 || size5 ~=0
                        if size4 > size5 
                            Sector2(i) = Sector4(i);
                        else
                            Sector2(i) = Sector5(i);
                        end
                 end
                 if Sector2(i) == 0
                        if size6 ~=0 || size8 ~=0
                         
                            if size6 > size8
                                Sector2(i) = Sector6(i);
                            else
                                Sector2(i) = Sector8(i);
                            end
                        else
                          Sector2(i) = Sector7(i);
                        end
                 end
             end
      end
       if Sector3(i) == 0
           if size2 ~=0 || size5 ~=0
             if size2 > size5
                 Sector3(i) = Sector2(i);
             else
                 Sector3(i) = Sector5(i);
             end
           end
           if Sector3(i) == 0
                if size1 ~=0 || size8 ~=0
                    if size1 > size8
                        Sector3(i) = Sector1(i);
                    else
                        Sector3(i) = Sector8(i);
                    end
                end
           end
                if Sector3(i) == 0
                    if size4 ~=0 || size7 ~=0
                        if size4 > size7 
                            Sector3(i) = Sector4(i);
                        else
                            Sector3(i) = Sector7(i);
                        end
                    else
                       Sector3(i) = Sector6(i);
                    end
                end
      end
       if Sector4(i) == 0
           if size6 ~=0 || size1 ~=0
             if size6 > size1
                 Sector4(i) = Sector6(i);
             else
                 Sector4(i) = Sector1(i);
             end
           end
             if Sector4(i) == 0
                 if size2 ~=0 || size7 ~=0
                      if size2 > size7 
                          Sector4(i) = Sector2(i);
                      else
                          Sector4(i) = Sector7(i);
                      end
                 end
             end
                if Sector4(i) == 0
                    if size3 ~=0 || size8 ~=0
                        if size3 > size8
                            Sector4(i) = Sector3(i);
                        else
                            Sector4(i) = Sector8(i);
                        end
                    else
                       Sector4(i) = Sector5(i);
                    end
                end
            
        end
      if Sector5(i) == 0
          if size3 ~=0 || size8 ~=0
             if size3 > size8
                 Sector5(i) = Sector3(i);
             else
                 Sector5(i) = Sector8(i);
             end
          end
             if Sector5(i) == 0
                 if size2 ~=0 || size7 ~=0
                      if size2 > size7
                         Sector5(i) = Sector2(i);
                      else
                         Sector5(i) = Sector7(i);
                      end
                 end
             end
                if Sector5(i) == 0
                    if size1 ~=0 || size6 ~=0
                         if size1 > size6
                             Sector5(i) = Sector1(i);
                         else
                             Sector5(i) = Sector6(i);
                         end
                     else
                       Sector5(i) = Sector4(i);
                    end
                end
      end
       if Sector6(i) == 0
           if size4 ~=0 || size7 ~=0
             if size4 > size7
                 Sector6(i) = Sector4(i);
             else
                 Sector6(i) = Sector7(i);
             end
           end
             if Sector6(i) == 0
                 if size1 ~=0 || size8 ~=0
                        if size1 > size8
                            Sector6(i) = Sector1(i);
                        else
                            Sector6(i) = Sector8(i);
                        end
                 end
             end
                if Sector6(i) == 0
                    if size2 ~=0 || size5 ~=0
                        if size2 > size5
                            Sector6(i) = Sector2(i);
                        else
                            Sector6(i) = Sector5(i);
                        end
                    else
                       Sector6(i) = Sector3(i);
                    end
                end
            
       end
       if Sector7(i) == 0
          if size6 ~=0 || size8 ~=0
             if size6 > size8
                 Sector7(i) = Sector6(i);
             else
                 Sector7(i) = Sector8(i);
             end
          end
             if Sector7(i) == 0
                 if size4 ~=0 || size5 ~=0
                     if size4 > size5
                         Sector7(i) = Sector4(i);
                     else
                         Sector7(i) = Sector5(i);
                     end
                 end
             end
                if Sector7(i) == 0
                    if size1 ~= 0 || size3 ~=0
                        if size1 > size3
                             Sector7(i) = Sector1(i);
                        else
                            Sector7(i) = Sector3(i);
                        end
                    else
                        Sector7(i) = Sector2(i);
                    end
                end
            
      end
      if Sector8(i) == 0
          if size5 ~=0 || size7 ~=0
             if size5 > size7
                 Sector8(i) = Sector5(i);
             else
                 Sector8(i) = Sector7(i);
             end
          end
             if Sector8(i) == 0
                 if size3 ~=0 || size6 ~=0
                      if size3 > size6
                         Sector8(i) = Sector3(i);
                      else
                         Sector8(i) = Sector6(i);
                      end
                 end
             end
                if Sector8(i) == 0
                   if size2 ~=0 || size4 ~=0
                        if size2 > size4
                           Sector8(i) = Sector2(i);
                        else
                           Sector8(i) = Sector4(i);
                        end
                   else
                       Sector8(i) = Sector1(i);
                   end
                end
      end
              
end
%% Final RA File Creation

for q=1:numRAgrains
    FinalRAFile(q,1) = grainfile(RAgrainLocation(q),1);  %ID
    FinalRAFile(q,2) = grainfile(RAgrainLocation(q),9);  %Area
    FinalRAFile(q,3) = grainfile(RAgrainLocation(q),10); %Diameter
    FinalRAFile(q,4) = grainfile(RAgrainLocation(q),11); %GSAR
    FinalRAFile(q,5) = grainfile(RAgrainLocation(q),12); %MAO
    FinalRAFile(q,6) = RASchmidFactor(q);                %RA Schmid Factor
    FinalRAFile(q,7) = grainfile(RAgrainLocation(q),13); %NumNeighbor
    FinalRAFile(q,8) = Sector1(q);
    FinalRAFile(q,9) = Sector2(q);
    FinalRAFile(q,10) = Sector3(q);
    FinalRAFile(q,11) = Sector4(q);
    FinalRAFile(q,12) = Sector5(q);
    FinalRAFile(q,13) = Sector6(q);
    FinalRAFile(q,14) = Sector7(q);
    FinalRAFile(q,15) = Sector8(q);
    FinalRAFile(q,16) = SurroundingFerritePercent(q);
    FinalRAFile(q,17) = SurroundingMartPercent(q);
    FinalRAFile(q,18) = SurroundingRApercent(q);
    
end

Headers = {'GrainID', 'Area', 'Diameter', 'GSAR', 'MAO', 'Schmid_Factor', 'Num_Neighbors','Sector_1','Sector_2','Sector_3','Sector_4','Sector_5','Sector_6','Sector_7','Sector_8', 'Surrounding_Ferrite', 'Surrounding_Martensite','Surrounding_RA'};

T = array2table(FinalRAFile,'VariableNames',Headers);
toDelete = T.Num_Neighbors<1; %Delete all rows (or RA grains) that have 0 neighbors, meaning they're on the edge of the scan or
T(toDelete,:) = [];

%Create CSV file for RA data by prompting user for filename. File will be
%saved in RA Grain Analysis folder regardless of the path created when
%making the filename.
savequery = menu('Save to a CSV File?', 'Yes','No');
if (savequery == 1)
   [filename,path] = uiputfile('*.csv','Save File As');
    writetable(T, filename);
end
%% Clear Useless Variables at end of file
clear cnames displaymax displaymin displaysuggestions f filename Headers i index j k m minlegend  neighborCols neighborRows numRARows numRACols numRAgrains
clear  path q r results sizeCols sizeRows t uilegend x row col n p

for i=1:29
    for j = 1:6
        ShieldingFactor(i,j) = NeighborPercentSurrounding(i,j)*RANeighborPhaseSchmid(i,j);
    end
end
    

