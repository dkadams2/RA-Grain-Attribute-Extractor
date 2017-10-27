%This code takes an exported map of values for ferrite and austenite Taylor
%Factors and adds them onto the grain file type 2. It locates the grains in
%the grain file and gives them the TF found at the same x,y coordinates in
%the map value txt file.
%Derrik Adams 10/26/14

clear all;

%Load Grain File
grainfilename = uigetfile({'*.txt'},'Select Grain File');
grainfile = dlmread(grainfilename);
[grainsizeRows, grainsizeCols] = size(grainfile);

%Load TF Map File
MapFilename = uigetfile({'*.txt'},'Select Map File');
mapfile = dlmread(MapFilename);
[mapsizeRows, mapsizeCols] = size(mapfile);

%Loop through the grain file to get the coordinates needed and then find
%those coordinates in the map file and assign the TF vector as such
for i=1:grainsizeRows
    Xpos = grainfile(i,5); %x-position of the grain 
    Ypos = grainfile(i,6); %y-position of the grain
    for j = 1:mapsizeRows
        
        %Need to make ranges on the variation of Xpos and Ypos of 1/2 the
        %step size in um since the map coordinates are by each point (so by
        %the step size) but the grain file coordinates are not on the same
        %step size, they are custom to position but fall in the correct
        %range...
        if  (mapfile(j,1)-.04) <= Xpos && Xpos <= (mapfile(j,1)+.04) && (mapfile(j,2)-.04) <= Ypos && Ypos <= (mapfile(j,2)+.04) 
            if mapfile(j,3) ~= -1
                TFgrainfile(i) = mapfile(j,3);
            else
                TFgrainfile(i) = mapfile(j,4);
            end
        end
    end
end


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