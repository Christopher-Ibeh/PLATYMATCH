%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            PLATY PARTICLE MATCHING CODE
%%                        BY
%%                   CHRISTOPHER IBEH
%%    DEPARTMENT OF CIVIL AND ENVIRONMENTAL ENGINEERING
%%              UNIVERSITY OF STRATHCLYDE
%%                      2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
%Import scan 1 and scan 2 particle sattributes from spreadsheet
%Spreadsheet columns and attributes stored:1- Particle index number;	
%2- Area;	3- Major axis length ;	4- Particle minor vertical axis orientation; 
%5- Particle minor horizontal axis orientation; 6- Particle major vertical axis orientation; 
%7- Particle major horizontal axis orientation	8- Particle intermediate axis length;	
%9- BaryCenterX;	10- BaryCenterY; 11- BaryCenterZ;	12- Perimeter	

%Import initial scan data from excel
data0=xlsread('DATA.xls','Scan_1'); 

%Import second scan data from excel
data1=xlsread('DATA.xls','Scan_2');

%% Prepare data for particle matching
%Check input data matrix size           
[m,n]=size(data0);
[p,o]=size(data1);

%count the number of data points
Q0 = numel(data0(:,1));
Q1 = numel(data1(:,1));

%Generate the particle numbers (i.e numbers from 1 to M in increaments of 1)
q0 = [1:Q0]';
q1 = [1:Q1];

%% Calculate normalised errors for ExtentMax1 (major principal axis), ExtentMax1 (Intermediate axis) , Perimeter and the sum of all the errors
%first define the search cylinder in scan 2 and around scan 1 particle centre

Z_conditionA=0 ; % downward displacement positive
Z_conditionB=-1; % upward displacement negative
radius_condition=0.4;
alldata=[];
  for j=1:m
     for i=1:p
              if (sqrt((data0(j,9) - data1(i,9)).^2 + (data0(j,10) - data1(i,10)).^2) <=radius_condition) &&...
                (((data0(j,11)-data1(i,11))>= Z_conditionB && (data0(j,11)-data1(i,11)) < 0) || ((data0(j,11)-data1(i,11)) <= Z_conditionA && (data0(j,11)-data1(i,11)) > 0))  &&...
                (((abs((data0(j,3)) - (data1(i,3)))/data0(j,3))+ (abs((data0(j,8)) - (data1(i,8)))/data0(j,8)) + (abs((data0(j,12)) - (data1(i,12)))/data0(j,12)))/3) <=0.05;
                
                alldata=[alldata; data0(j,1:12), data1(i,1:12)];
    
              else continue;
              
              end  
      end
 end
  
%% Calculate total differences of all measures as well as displacements
% Will be used to determine correct match
[r,s]=size(alldata);

 for j=1:r
    alldata(j,25)= sqrt((alldata(j,9)-alldata(j,21))^2+(alldata(j,10)-alldata(j,22))^2 +(alldata(j,11)-alldata(j,23))^2);% Eucledean Displacement
    alldata(j,26)=(alldata(j,21)-alldata(j,9)); %displacement in the x direction
    alldata(j,27)=(alldata(j,22)-alldata(j,10)); %displacement in the y direction
    alldata(j,28)=(alldata(j,23)-alldata(j,11)); %displacement in the z direction    
%rotation
    alldata(j,29)=(alldata(j,16)-alldata(j,4)); %minor axis rotation
    alldata(j,30)=(alldata(j,18)-alldata(j,6)); %major axis rotation
    alldata(j,31)=(((abs((alldata(j,3)) - (alldata(j,15)))/alldata(j,3))+ (abs((alldata(j,8)) - (alldata(j,20)))/alldata(j,8)) + (abs((alldata(j,12)) - (alldata(j,24)))/alldata(j,12)))/3);
 end
 
%% Remove duplicates (more than 1 match) by choosing particle with least matching error value
% Check both the first particle and second Particle

[y,z] = size(alldata);

tic
for i=1:y
    ind=find(alldata(:,1)==alldata(i,1));
    if alldata(i,1)==missing; continue;
    elseif length(ind)==1; continue;
    else

        ind2=find(alldata(ind,31)==min(alldata(ind,31)));
        ind(ind2)=[];
        alldata(ind,:)=NaN;
        
    end
    ind=[]; ind2=[];
end

alldata(any(isnan(alldata),2),:)=[];

[y,z] = size(alldata);
for i=1:y %length(alldata)
    ind=find(alldata(:,13)==alldata(i,13));
    if alldata(i,19)==missing; continue;
    elseif length(ind)==1; continue;
    else

        ind2=find(alldata(ind,31)==min(alldata(ind,31)));
        ind(ind2)=[];
        alldata(ind,:)=NaN;
        
    end
    ind=[]; ind2=[];
end
toc
alldata(any(isnan(alldata),2),:)=[];
 
%% Prepare data for plots
%Extract change in x-displacement and store in DX
DX=alldata(:,26);
%Extract change in y-displacement and store in DY
DY=alldata(:,27);
%Extract change in x-displacement and store in DZ
DZ=alldata(:,28);

%Store the initial and Finalparticle positions X,Y,Z coordinates as BCx1,BCy1,BCz1 and BCx1,BCy1,BCz1,also store the eucledian distance for each particle between scans as w. 
BCx1=alldata(:,9);
BCy1=alldata(:,10);
BCz1=alldata(:,11);
BCx2=alldata(:,21);
BCy2=alldata(:,22);
BCz2=alldata(:,23);
w=alldata(:,25); 

%Store particle rotation data
u=alldata(:,29);% Minor axis rotation
v=alldata(:,30);% Major axis rotation
s = ones(size(alldata,1), 1) * 2; %sizes of markers

%% Plot figures
figure(1) %Plot Euclidean distance map
scatter3(BCx1,BCy1,BCz1,s,w,'filled');
colorbar;
xlabel('Particle x Position (mm)');
ylabel('Particle y Position (mm)');
zlabel('Particle z Position (mm)');
title('Euclidean distance (mm)');

figure(2) %Plot z- displacement map 
scatter3(BCx1,BCy1,BCz1,s,DZ,'filled');
xlabel('Particle y Position (mm)');
ylabel('Particle z Position (mm)');
zlabel('Particle z Position (mm)');
colorbar;
title('Displacement in Z direction (mm)');
xlim([0 12])
ylim([0 12])
zlim([0 25])

figure(3)%Plot horizontal rotation map
scatter3(BCx1,BCy1,BCz1,s,u,'filled');
colorbar;
xlabel('Particle y Position (mm)');
ylabel('Particle z Position (mm)');
zlabel('Particle z Position (mm)');
title('ROTATION');
xlim([0 12])
ylim([0 12])
zlim([0 25])

figure(4)%Plot displacement vector map
quiver3(BCx1,BCy1,BCz1,DX,DY,DZ);
xlabel('Particle y Position (mm)');
ylabel('Particle z Position (mm)');
zlabel('Particle z Position (mm)');
title('Map of Particle displacement vector');
xlim([0 12])
ylim([0 12])
zlim([0 25])

figure(6)% Plot Particle z position Vs Particle z-displacement 
[p,S] = polyfit(BCz1,DZ,1); 
[DZ_fit,delta] = polyval(p,BCz1,S);
plot(BCz1,DZ,'bo')
hold on
plot(BCz1,DZ_fit,'r-')
legend('Data','- Average Displacement')
xlabel('Particle z Position');
ylabel('Particle z-displacement (mm)');
title('Particle z position Vs particle z-displacement'); 
grid on;
% Place equation of the average displacement in upper left of graph.
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.05 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('y = %f * x + %f', p(1), p(2));
text(xt, yt, caption, 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');