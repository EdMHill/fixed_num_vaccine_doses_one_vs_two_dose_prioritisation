%% Script to investigate prioritisation of a one dose or two dose vaccination 
%% schedule given a fixed number of vaccine doses and with respect 
%% to a measure of maximising averted deaths

clear

%% LOAD POPULATION SIZES & RELATIVE RISKS FOR EACH PRIORITY GROUP
load group_popns_and_rel_risks.mat
   % RPP - Population in each group (Priority Group estimate)
   % RPPa - Population in each group (Age Only estimate)
   % RRisk - Per group, deaths per person (Priority Group estimate)
   % RRiska - Per group, deaths per person (Age Only estimate)

%% HOMOGENOUS ALLOCATION STRATEGY

% Proportion of each age group that gets vaccinated. 
p=0.9;     

PP1= 0.45:0.001:0.99;  % Relative strength of first dose compared to having 2 doses.
DD = 1e4:1e4:36e6;      % Number of doses to examine.

% Initialise array to compare reduction in risks of one dose priority
% vs two dose priority
Ratio = zeros(length(PP1),length(DD));

%%% Perform calculation for "Priority Group Estimate" %%%
% Iterate over each relative first dose efficacy & amount of doses
for x=1:length(PP1)
    p1=PP1(x);
    
    for y=1:length(DD)
        D=DD(y);
        
        % Initialise variables tracking number of vaccines allocated
        % in one dose priority strategy
        % Vectors: Entry per age-group.
        V1=0*RPP;
        V12=V1;   % Second doses administed in first dose prioritised strategy
        
        % Initialise variables tracking number of vaccines allocated
        % in two dose priority strategy
        % Vector: Entry per age-group.
        V2=V1;   
        
        % One dose priority loop
        d=D;
        for i=1:9
            % Iterate over each priority group
            % Check if doses remain to be allocated
            if RPP(i)*p<d
                % Doses allocated to group i. Reduce number of doses that remain to be allocated.
                V1(i)=RPP(i)*p;
                d=d-RPP(i)*p;   
            else
                % All remaining doses/no doses allocated to group i
                V1(i)=d;
                d=0;
            end
        end
        
        % All groups checked for first dose assignment.
        % If any doses remain, can be allocated as second doses.
        if d>0
            for i=1:9
                if RPP(i)*p<d
                    % Doses allocated to group i. Reduce number of doses that remain to be allocated.
                    V12(i)=RPP(i)*p;
                    d=d-RPP(i)*p;
                else
                    % All remaining doses/no doses allocated to group i
                    V12(i)=d;
                    d=0;
                end
            end
        end
        
        % Compute reduction in risk
        % First term: Two dose efficacy; Second term: One dose efficacy.
        P1=RRisk.*V12 + RRisk.*(V1-V12).*p1;
        
        % Two dose priority strategy
        d=D/2;
        for i=1:9
            if RPP(i)*p<d
                V2(i)=RPP(i)*p;
                d=d-RPP(i)*p;
            else
                V2(i)=d;
                d=0;
            end
        end
        
        % Compute reduction in risk. All two dose efficacy.
        P2=RRisk.*V2;
        
        % Allocate ratio of risk reduction to storage array
        Ratio(x,y)=sum(P1)/sum(P2);
    end
end

%%% Perform calculation using age-only estimate %%%

clear Ratioa

% Initialise array to compare reduction in risks of one dose priority
% vs two dose priority
Ratioa = zeros(length(PP1),length(DD));

% Iterate over each relative first dose efficacy & amount of doses
for x=1:length(PP1)
    p1=PP1(x);
    
    for y=1:length(DD)
        D=DD(y);
                
        % Initialise variables tracking number of vaccines allocated
        % in one dose priority strategy
        % Vectors: Entry per age-group.
        V1=0*RPPa;
        V12=V1;   % Second doses administed in first dose prioritised strategy
        
        % Initialise variables tracking number of vaccines allocated
        % in two dose priority strategy
        % Vector: Entry per age-group.
        V2=V1; 
        
        % One dose priority loop
        d=D;
        for i=1:9
             % Iterate over each priority group
            % Check if doses remain to be allocated
            if RPPa(i)*p<d
                % Doses allocated to group i. Reduce number of doses that remain to be allocated.
                V1(i)=RPPa(i)*p;
                d=d-RPPa(i)*p;
            else
                % All remaining doses/no doses allocated to group i
                V1(i)=d;
                d=0;
            end
        end
        
        % All groups checked for first dose assignment.
        % If any doses remain, can be allocated as second doses.
        if d>0
            for i=1:9
                if RPPa(i)*p<d
                    % Doses allocated to group i. Reduce number of doses that remain to be allocated.
                    V12(i)=RPPa(i)*p;
                    d=d-RPPa(i)*p;
                else
                    % All remaining doses/no doses allocated to group i
                    V12(i)=d;
                    d=0;
                end
            end
        end
        
        % Compute reduction in risk
        % First term: Two dose efficacy; Second term: One dose efficacy.
        P1=RRiska.*V12 + RRiska.*(V1-V12).*p1;
        
        % Two dose priority strategy
        d=D/2;
        for i=1:9
            if RPPa(i)*p<d
                % Doses allocated to group i. Reduce number of doses that remain to be allocated.
                V2(i)=RPPa(i)*p;
                d=d-RPPa(i)*p;
            else
                % All remaining doses/no doses allocated to group i
                V2(i)=d;
                d=0;
            end
        end
        
        % Compute reduction in risk. All two dose efficacy.
        P2=RRiska.*V2;
        
        % Allocate ratio of risk reduction to storage array
        Ratioa(x,y)=sum(P1)/sum(P2);
    end
end

%% PLOT THE RESULTS FOR HOMOGENEOUS STRATEGY
close all

% Priority Group Estimate
figure(1);
colormap([0.9 0.9 0.9 ; 0.7 0.7 0.7]);
pcolor(DD/1e6,100*PP1,sign(Ratio-1.0001)); shading flat;
hold on
%plot(p*sum(RPP)/1e6+[0 0],[0 100],'k');
%[c,h]=contour(DD/1e6,100*PP1,Ratioa,[1 1]);
%set(h,'LineWidth',2,'Color','c');
[~,h]=contour(DD/1e6,100*PP1,Ratio,[1 1]);
set(h,'LineWidth',2,'Color','m');
hold off
text(15,60,'Prioritise Second Dose','HorizontalAlignment','center','FontSize',16);
text(6,95,'Prioritise First Dose','HorizontalAlignment','left','FontSize',16);
xlabel('Total Vaccine Doses (millions)','FontSize',14);
ylabel({'Relative Efficacy from First Dose(%)',' '},'FontSize',14);
axis([0 p*sum(RPP)/1e6 45 99]);

set(gca,'FontSize',16)
set(gca,'LineWidth',1)
box on

print -dpng Homogeneous_Dose_Priority.png

% Age only estimates
figure(2);
colormap([0.9 0.9 0.9 ; 0.7 0.7 0.7]);
pcolor(DD/1e6,100*PP1,sign(Ratioa-1.0001)); shading flat;
hold on
%plot(p*sum(RPPa)/1e6+[0 0],[0 100],'k');
[~,h]=contour(DD/1e6,100*PP1,Ratioa,[1 1]);
set(h,'LineWidth',2,'Color','c');
%[c,h]=contour(DD/1e6,100*PP1,Ratio,[1 1]);
%set(h,'LineWidth',2,'Color','m');
hold off
text(13,70,'Prioritise Second Dose','HorizontalAlignment','center','FontSize',16);
text(4,95,'Prioritise First Dose','HorizontalAlignment','left','FontSize',16);
xlabel('Total Vaccine Doses (millions)','FontSize',14);
ylabel({'Relative Efficacy from First Dose(%)',' '},'FontSize',14);
axis([0 p*sum(RPPa)/1e6 45 99]);

set(gca,'FontSize',16)
set(gca,'LineWidth',1)
box on

print -dpng Homogeneous_Dose_Priority_AgeOnly.png

%% HETEROGENEOUS ALLOCATION STRATEGY: AGE ONLY ESTIMATE
% Set relative first dose efficacies to be assessed
relVE = 0.45:0.01:0.99;

% Set amount of doses to be assessed
TV = [0:0.01:24]*1e6; 
clear TV1 TV2   % TV1 & TV2: Total vaccine doses used as first dose and second dose, respectively.

% Uptake proportion (p) & absolute numbers who receive vaccine doses per
% age group
p=0.9;
pRPP=p*RPPa;

% Iterate over each vaccine efficacy
% For given number of doses, determine split between first and second dose
% allocation
for rv=1:length(relVE)
    
    % Initialise tracking variables
    % Index of groups that have been covered by first dose (v1), second
    % dose (v2)
    v1=0; v2=0; 
    
    % Intialise tracking variables, storing amount of vaccine doses
    % allocated after each iteration of allocation loop
    PX=[0]; V1=[0]; V2=[0];
    
    % Run allocation loop
    for z=1:10
        
        % Increment group index to work out which next step is best.
        % t1 - Index for the next group to receive first dose
        % t2 - Index for the next group to receive second dose
        t1=v1+1;   t2=v2+1;
        
        % Check selected groups have population eligible for uptake
        % If not, increment group index (t1,t2) by one.
        if t1<10 
            if pRPP(t1)==0  
                t1=t1+1; 
            end
        end
        
        if pRPP(t2)==0  
            t2=t2+1; 
        end
        
        % Test adding v2 (if seond dose allocation behind first dose
        % allocation in terms of groups covered)
        if v2<v1
            % Use 1-relVE(rv), the added efficacy from second dose
            % conditional on having had first dose
            Benefit2=RRiska(t2)*(1-relVE(rv));
        else
            Benefit2=-10;
        end
        
        % Check if all groups have received first dose.
        if t1>9
            Benefit1=0;
        else % If not, compute benefit
            Benefit1=RRiska(t1)*relVE(rv);
        end
        
        % Check benefits between two options
        if Benefit1 > Benefit2 % best to give a first dose to group t1
            % [Vaccine doses previously allocated  Updated amount of vaccine dose previously allocated]
            X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:t1))+sum(pRPP(1:v2))]/1e6;
            
            % Append doses allocated to tracking vectors
            PX=[PX sum(pRPP(1:t1))+sum(pRPP(1:v2))]; V1=[V1 sum(pRPP(1:t1))]; V2=[V2 sum(pRPP(1:v2))];
            
            % Update v1, tracking groups covered by first doses
            v1=t1; 
        else % best to give a second dose to group t2
            % [Vaccine doses previously allocated  Updated amount of vaccine dose previously allocated]
            X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:v1))+sum(pRPP(1:t2))]/1e6;
            
            % Append doses allocated to tracking vectors
            PX=[PX sum(pRPP(1:v1))+sum(pRPP(1:t2))]; V1=[V1 sum(pRPP(1:v1))]; V2=[V2 sum(pRPP(1:t2))];
            
            % Update v2, tracking groups covered by second dose
            v2=t2; 
        end
        
        % Return number of vaccines allocated (split by first dose & second dose) 
        % at total dose availability values of interest
        TV1(:,rv)=interp1(PX,V1,TV);
        TV2(:,rv)=interp1(PX,V2,TV);
    end
end

% Manually set for the scenarios where relative efficacy of first dose is
% below 50%, that number of doses is equally split between first and second
% doses
TV1(:,relVE<0.5)=TV'*ones(1,sum(relVE<0.5))/2;
TV2(:,relVE<0.5)=TV'*ones(1,sum(relVE<0.5))/2;

TV1(isnan(TV1))=1;  TV2(isnan(TV2))=1;
  
%% PLOT HEATMAP HETEROGENEOUS ALLOCATION STRATEGY: AGE ONLY ESTIMATE
close all
figure(1);

pcolor(TV/1e6,100*relVE,100*TV2'./(TV1'+TV2')); shading flat; caxis([0 50]); colorbar;
hold on
plot([0 100],[70 70],'-k','LineWidth',1,'Color',[0.99 0.99 0.99]);
plot([0 100],[80 80],'-k','LineWidth',1,'Color',[0.99 0.99 0.99]);
plot([0 100],[90 90],'-k','LineWidth',1,'Color',[0.99 0.99 0.99]);
plot([0 100],[70 70],'-k');
plot([0 100],[80 80],'-k');
plot([0 100],[90 90],'-k');

Ratioa(:,DD>sum(p*RPPa))=NaN;

[c,h]=contour(DD/1e6,100*PP1,Ratioa,[1 1]);
set(h,'LineWidth',3,'Color','k');
[c,h]=contour(DD/1e6,100*PP1,Ratioa,[1 1]);
set(h,'LineWidth',2,'Color','c');

plot(sum(p*RPPa)/1e6+[0 0],[0 100],'--k','LineWidth',2);
plot(sum(p*RPPa)/1e6+[0 0],[0 100],'--c','LineWidth',1);

hold off

xlabel('Total Vaccine Doses (millions)','FontSize',14);
ylabel({'Relative Efficacy from First Dose(%)',' '},'FontSize',14);
H=colorbar;
set(get(H,'YLabel'),'String','Optimal proportion of vaccines used for second dose','FontSize',16);

axis([0 24 45 99]);
set(gca,'Fontsize',16)
box on

print -dpng Heterogeneous_Dose_Priority_AgeOnly.png

%% HETEROGENEOUS ALLOCATION STRATEGY: PRIORITY GROUP ESTIMATE

% Set relative first dose efficacies to be assessed
relVE=[0.45:0.01:0.99];

% Set amount of doses to be assessed
TV=[0:0.01:36]*1e6; clear TV1 TV2

% Uptake proportion (p) & absolute numbers who receive vaccine doses per
% age group
p=0.9;
pRPP=p*RPP;

% Iterate over each vaccine efficacy
% For given number of doses, determine split between first and second dose
% allocation
for rv=1:length(relVE)
    
    % Initialise tracking variables
    % Index of groups that have been covered by first dose (v1), second
    % dose (v2)
    v1=0; v2=0; 
    
    % Intialise tracking variables, storing amount of vaccine doses
    % allocated after each iteration of allocation loop
    PX=[0]; V1=[0]; V2=[0];
    
    % Run allocation loop
    for z=1:16
        
        % Increment group index to work out which next step is best.
        % t1 - Index for the next group to receive first dose
        % t2 - Index for the next group to receive second dose
        t1=v1+1;   t2=v2+1;
        
        % Check selected groups have population eligible for uptake
        % If not, increment group index (t1,t2) by one.
        if t1<10 
            if pRPP(t1)==0  
                t1=t1+1; 
            end
        end
        
        if pRPP(t2)==0  
            t2=t2+1; 
        end
        
        % Test adding v2 (if seond dose allocation behind first dose
        % allocation in terms of groups covered)
        if v2<v1
            % Use 1-relVE(rv), the added efficacy from second dose
            % conditional on having had first dose
            Benefit2=RRisk(t2)*(1-relVE(rv));
        else
            Benefit2=-10;
        end
        
        % Check if all groups have received first dose.
        if t1>9
            Benefit1=0;
        else % If not, compute benefit
            Benefit1=RRisk(t1)*relVE(rv);
        end
        
        % Check benefits between two options
        if Benefit1 > Benefit2 % best to give a first dose to group t1
            % [Vaccine doses previously allocated  Updated amount of vaccine dose previously allocated]
            X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:t1))+sum(pRPP(1:v2))]/1e6;
            
            % Append doses allocated to tracking vectors
            PX=[PX sum(pRPP(1:t1))+sum(pRPP(1:v2))]; V1=[V1 sum(pRPP(1:t1))]; V2=[V2 sum(pRPP(1:v2))];
            
            % Update v1, tracking groups covered by first doses
            v1=t1; 
        else % best to give a second dose to group t2
            % [Vaccine doses previously allocated  Updated amount of vaccine dose previously allocated]
            X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:v1))+sum(pRPP(1:t2))]/1e6;
            
            % Append doses allocated to tracking vectors
            PX=[PX sum(pRPP(1:v1))+sum(pRPP(1:t2))]; V1=[V1 sum(pRPP(1:v1))]; V2=[V2 sum(pRPP(1:t2))];
            
            % Update v2, tracking groups covered by second dose
            v2=t2; 
        end
        
        % Return number of vaccines allocated (split by first dose & second dose) 
        % at total dose availability values of interest
        TV1(:,rv)=interp1(PX,V1,TV);
        TV2(:,rv)=interp1(PX,V2,TV);
    end
end


% Manually set for the scenarios where relative efficacy of first dose is
% below 50%, that number of doses is equally split between first and second
% doses
TV1(:,relVE<0.5)=TV'*ones(1,sum(relVE<0.5))/2;
TV2(:,relVE<0.5)=TV'*ones(1,sum(relVE<0.5))/2;
  
TV1(isnan(TV1))=1;  TV2(isnan(TV2))=1;
  
%% PLOT HEATMAP HETEROGENEOUS ALLOCATION STRATEGY: PRIORITY GROUP ESTIMATE
close all
figure(1);

pcolor(TV/1e6,100*relVE,100*TV2'./(TV1'+TV2')); shading flat; caxis([0 50]); colorbar;
hold on
plot([0 100],[70 70],'-k','LineWidth',1,'Color',[0.99 0.99 0.99]);
plot([0 100],[80 80],'-k','LineWidth',1,'Color',[0.99 0.99 0.99]);
plot([0 100],[90 90],'-k','LineWidth',1,'Color',[0.99 0.99 0.99]);
plot([0 100],[70 70],'-k');
plot([0 100],[80 80],'-k');
plot([0 100],[90 90],'-k');


Ratio(:,DD>sum(p*RPP))=NaN;

[c,h]=contour(DD/1e6,100*PP1,Ratio,[1 1]);
set(h,'LineWidth',3,'Color','k');
[c,h]=contour(DD/1e6,100*PP1,Ratio,[1 1]);
set(h,'LineWidth',2,'Color','m');

plot(sum(p*RPP)/1e6+[0 0],[0 100],'--k','LineWidth',2);
plot(sum(p*RPP)/1e6+[0 0],[0 100],'--m','LineWidth',1);

hold off

xlabel('Total Vaccine Doses (millions)','FontSize',14);
ylabel({'Relative Efficacy from First Dose(%)',' '},'FontSize',14);
H=colorbar;
set(get(H,'YLabel'),'String','Optimal proportion of vaccines used for second dose','FontSize',16);

axis([0 36 45 99]);
set(gca,'Fontsize',16)
box on

print -dpng Heterogeneous_Dose_Priority.png

%% DOSES OVER TIME (Age only estimate)

% Plot colour settings
C=[0 113 187; 217 83 25; 237 177 32; 126 47 142; 119 172 48; 77 190 238; 162 20 47; 100 100 100; 180 255 180 ]/255;
TC=[1 1 1; 0 0 0; 0 0 0; 1 1 1; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 0 0 0];

% Uptake proportion (p) & absolute numbers who receive vaccine doses per
% age group
p=0.9;
pRPP=p*RPPa;

% Set the relative vaccine efficacies to be tested
relVE=[0.7 0.8 0.9];
clf; set(gcf,'position',[45 743 890 602]);

MX=24; % Maximum number of doses being considered

for rv=1:length(relVE)
    
    subplot(length(relVE),1,rv);
    
    % Perform heterogeneous allocation strategy
    % Add next group selected to plot after each iteration of loop
    
    % Initialise tracking variables
    % Index of groups that have been covered by first dose (v1), second
    % dose (v2)
    v1=0; v2=0;
    
    % Run allocation loop
    for z=1:10
        
        % Increment group index to work out which next step is best.
        % t1 - Index for the next group to receive first dose
        % t2 - Index for the next group to receive second dose
        t1=v1+1;   t2=v2+1;
        
        % Check selected groups have population eligible for uptake
        % If not, increment group index (t1,t2) by one.
        if t1<10 
            if pRPP(t1)==0  
                t1=t1+1; 
            end 
        end
        
        if pRPP(t2)==0  
            t2=t2+1; 
        end
        
        % test adding v2
        if v2<v1 
            Benefit2=RRiska(t2)*(1-relVE(rv));
        else
            Benefit2=-10;
        end
        
        % Check if all groups have received first dose.
        if t1>9
            Benefit1=0;
        else
            Benefit1=RRiska(t1)*relVE(rv);
        end
        
       % Check benefits between two options
        if Benefit1 > Benefit2 % best to give a first dose to t1 !
            % [Vaccine doses previously allocated  Updated amount of vaccine dose previously allocated]
            X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:t1))+sum(pRPP(1:v2))]/1e6;
            fill(X([1 2 2 1]),1+0.3*[-1 -1 1 1],'r','FaceColor',C(t1,:)); hold on
            if mean(X)<MX
                text(mean(X),1,num2str(t1),'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','Color',TC(t1,:));
            end
            
            % Update v1, tracking groups covered by first doses
            v1=t1;
        else
            % [Vaccine doses previously allocated  Updated amount of vaccine dose previously allocated]
            X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:v1))+sum(pRPP(1:t2))]/1e6;
            fill(X([1 2 2 1]),2+0.3*[-1 -1 1 1],'r','FaceColor',C(t2,:)); hold on
            if mean(X)<MX
                text(mean(X),2,num2str(t2),'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','Color',TC(t2,:));
            end
            
            % Update v2, tracking groups covered by second dose
            v2=t2;
        end
        drawnow;
    
    end
    hold off
    
    % Plot properties
    axis([0 MX 0.5 2.5]);
    title(['Relative Efficacy ' num2str(100*relVE(rv)) '%']);
    
    if rv<length(relVE)
        set(gca,'XTick',[0:4:MX]);
    else
        xlabel({'Doses in millions'},'FontSize',16);
        set(gca,'XTick',[0:4:MX]);  
    end
    
   set(gca,'YTick',[1 2],'YTickLabel',{'Dose 1','Dose 2'},'FontSize',14);

end
    
print -dpng Doses_over_Time_AgeOnly.png

%% DOSES OVER TIME (Priority group estimate)

% Plot colour settings
C=[0 113 187; 217 83 25; 237 177 32; 126 47 142; 119 172 48; 77 190 238; 162 20 47; 100 100 100; 180 255 180 ]/255;
TC=[1 1 1; 0 0 0; 0 0 0; 1 1 1; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 0 0 0];

% Uptake proportion (p) & absolute numbers who receive vaccine doses per
% age group
p=0.9;
pRPP=p*RPP;

% Set the relative vaccine efficacies to be tested
relVE=[0.7 0.8 0.9];
clf; set(gcf,'position',[45 743 890 602]);
MX=36; % Maximum number of doses being considered

for rv=1:length(relVE)
    
    subplot(length(relVE),1,rv);
    
    % Perform heterogeneous allocation strategy
    % Add next group selected to plot after each iteration of loop
    
    % Initialise tracking variables
    % Index of groups that have been covered by first dose (v1), second
    % dose (v2)
    v1=0; v2=0;
    
    % Run allocation loop
    for z=1:16
        
        % Increment group index to work out which next step is best.
        % t1 - Index for the next group to receive first dose
        % t2 - Index for the next group to receive second dose
        t1=v1+1;   t2=v2+1;
        
        % Check selected groups have population eligible for uptake
        % If not, increment group index (t1,t2) by one.
        if t1<10 
            if pRPP(t1)==0  
                t1=t1+1; 
            end
        end
        if pRPP(t2)==0  
            t2=t2+1;
        end
        
        % Test adding v2 (if seond dose allocation behind first dose
        % allocation in terms of groups covered)
        if v2<v1 % test adding v2
            Benefit2=RRisk(t2)*(1-relVE(rv));
        else
            Benefit2=-10;
        end
        
        % Check if all groups have received first dose.
        if t1>9
            Benefit1=0;
        else
            Benefit1=RRisk(t1)*relVE(rv);
        end
        
        if Benefit1 > Benefit2 % best to give a first dose to t1 !
            % [Vaccine doses previously allocated  Updated amount of vaccine dose previously allocated]
            X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:t1))+sum(pRPP(1:v2))]/1e6;
            fill(X([1 2 2 1]),1+0.3*[-1 -1 1 1],'r','FaceColor',C(t1,:)); hold on
            if mean(X)<MX
                text(mean(X),1,num2str(t1),'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','Color',TC(t1,:));
            end
            
            % Update v1, tracking groups covered by second dose
            v1=t1;
        else
            % [Vaccine doses previously allocated  Updated amount of vaccine dose previously allocated]
            X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:v1))+sum(pRPP(1:t2))]/1e6;
            fill(X([1 2 2 1]),2+0.3*[-1 -1 1 1],'r','FaceColor',C(t2,:)); hold on
            if mean(X)<MX
                text(mean(X),2,num2str(t2),'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','Color',TC(t2,:));
            end
            
            % Update v2, tracking groups covered by second dose
            v2=t2;
        end
    drawnow;
    
    end
    hold off
    
    % Plot properties
    axis([0 MX 0.5 2.5]);
    title(['Relative Efficacy ' num2str(100*relVE(rv)) '%']);
    
    if rv<length(relVE)
        set(gca,'XTick',0:4:MX);
    else
        xlabel({'Doses in millions'},'FontSize',16);
        set(gca,'XTick',0:4:MX);
    end
    
    set(gca,'YTick',[1 2],'YTickLabel',{'Dose 1','Dose 2'},'FontSize',14);
    
end
    
print -dpng Doses_over_Time.png

%% BAR GRAPHS

% Plot colour settings
C=[0 113 187; 217 83 25; 237 177 32; 126 47 142; 119 172 48; 77 190 238; 162 20 47; 100 100 100; 180 255 180 ]/255;
TC=[1 1 1; 0 0 0; 0 0 0; 1 1 1; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 0 0 0];

% Uptake proportion (p) & absolute numbers who receive vaccine doses per
% age group
p=0.9;
pRPPa=p*RPPa;
pRPP=p*RPP;

% Set the relative vaccine efficacies to be tested
relVE=[0.7 0.8 0.9];
clf; set(gcf,'position',[64   300   827   983]);

for rv=1:length(relVE)
    
    % Specify the amount of doses (millions)
    DD=[4:4:24];
    
    % Iterate over the collection of overall dose amounts
    for D=1:length(DD)
        subplot(3,7,D+(rv-1)*7);
        
        %%% Priority group estimate %%%
        
        % Initialise tracking variables
        % Index of groups that have been covered by first dose (v1), second
        % dose (v2)
        v1=0; v2=0;
        
        % Intialise tracking variables, storing amount of vaccine doses
        % allocated after each iteration of allocation loop
        PX=[0]; V1=[0]; V2=[0];
        
        % Initialise variables to store max dose allocation amounts
        Mx1=0; Mx2=0; Mx1a=0; Mx2a=0;

        % Run allocation loop
        for z=1:16
            
            % Increment group index to work out which next step is best.
            % t1 - Index for the next group to receive first dose
            % t2 - Index for the next group to receive second dose
            t1=v1+1;   t2=v2+1;
            
            % Check selected groups have population eligible for uptake
            % If not, increment group index (t1,t2) by one.
            if t1<10 
                if pRPP(t1)==0  
                    t1=t1+1; 
                end
            end
            if pRPP(t2)==0  
                t2=t2+1; 
            end
            
            % Test adding v2 (if seond dose allocation behind first dose
            % allocation in terms of groups covered)
            if v2<v1
                Benefit2=RRisk(t2)*(1-relVE(rv));
            else
                Benefit2=-10;
            end
            
            % Check if all groups have received first dose.
            if t1>9
                Benefit1=0;
            else
                Benefit1=RRisk(t1)*relVE(rv);
            end
            
            % Check benefits between two options
            if Benefit1 > Benefit2 % best to give a first dose to t1 !
                X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:t1))+sum(pRPP(1:v2))]/1e6;
                PX=[PX sum(pRPP(1:t1))+sum(pRPP(1:v2))]; V1=[V1 sum(pRPP(1:t1))]; V2=[V2 sum(pRPP(1:v2))];
                if PX(end-1)<DD(D)*1e6  % below the limit so plot
                    if PX(end)<DD(D)*1e6  % won't exceed the limit
                        fill(2+0.4*[-1 1 1 -1],V1(end+[-1 -1 0 0])/1e6,'r','FaceColor',C(t1,:)); hold on; drawnow;
                        Mx1=V1(end)/1e6;
                    else % Case where does exceed the limit
                        tmp=V1; tmp(end)=tmp(end)-(PX(end)-DD(D)*1e6);
                        fill(2+0.4*[-1 1 1 -1], tmp(end+[-1 -1 0 0])/1e6 ,'r','FaceColor',C(t1,:)); hold on; drawnow;
                        Mx1=tmp(end)/1e6;
                    end
                end
                v1=t1;
            else % best to give a second dose to group t2
                X=[sum(pRPP(1:v1))+sum(pRPP(1:v2)) sum(pRPP(1:v1))+sum(pRPP(1:t2))]/1e6;
                PX=[PX sum(pRPP(1:v1))+sum(pRPP(1:t2))]; V1=[V1 sum(pRPP(1:v1))]; V2=[V2 sum(pRPP(1:t2))];
                if PX(end-1)<DD(D)*1e6  % below the limit so plot
                    if PX(end)<DD(D)*1e6  % won't exceded the limit
                        fill(5+0.4*[-1 1 1 -1],V2(end+[-1 -1 0 0])/1e6,'r','FaceColor',C(t2,:)); hold on; drawnow;
                        Mx2=V2(end)/1e6;
                    else % Case where does exceed the limit
                        tmp=V2; tmp(end)=tmp(end)-(PX(end)-DD(D)*1e6);
                        fill(5+0.4*[-1 1 1 -1], tmp(end+[-1 -1 0 0])/1e6,'r','FaceColor',C(t2,:)); hold on; drawnow;
                        Mx2=tmp(end)/1e6;
                    end
                end
                v2=t2;
            end
            
            % Plot properties
            if D>1
                set(gca,'YTickLabel',[]);
            else
                ylabel({'Optimal Deployment','millions of doses'},'FontSize',12);
                text(0,25,['Relative Vaccine Efficacy of First Dose = ' num2str(relVE(rv)*100) '%'],'FontWeight','Bold','FontSize',12);
            end
            
        end
        
        %%% Age only estimate %%%
        
        % Reinitialise variables
        v1=0; v2=0; PX=[0]; V1=[0]; V2=[0];

        % Run allocation loop
        for z=1:14
            
            % Increment group index to work out which next step is best.
            % t1 - Index for the next group to receive first dose
            % t2 - Index for the next group to receive second dose
            t1=v1+1;   t2=v2+1;
            
            % Check selected groups have population eligible for uptake
            % If not, increment group index (t1,t2) by one.
            if t1<10 if pRPPa(t1)==0  t1=t1+1; end; end
            if pRPPa(t2)==0  t2=t2+1; end
            
            % Test adding v2 (if seond dose allocation behind first dose
            % allocation in terms of groups covered)
            if v2<v1
                Benefit2=RRiska(t2)*(1-relVE(rv));
            else
                Benefit2=-10;
            end
            
            % Check if all groups have received first dose.
            if t1>9
                Benefit1=0;
            else
                Benefit1=RRiska(t1)*relVE(rv);
            end
            
            % Check benefits between two options
            if Benefit1 > Benefit2 % best to give a first dose to t1 !
                X=[sum(pRPPa(1:v1))+sum(pRPPa(1:v2)) sum(pRPPa(1:t1))+sum(pRPPa(1:v2))]/1e6;
                PX=[PX sum(pRPPa(1:t1))+sum(pRPPa(1:v2))]; V1=[V1 sum(pRPPa(1:t1))]; V2=[V2 sum(pRPPa(1:v2))];
                if PX(end-1)<DD(D)*1e6  % below the limit so plot
                    if PX(end)<DD(D)*1e6  % won't exceded the limit
                        fill(1+0.4*[-1 1 1 -1],V1(end+[-1 -1 0 0])/1e6,'r','FaceColor',C(t1,:)); hold on; drawnow;
                        Mx1a=V1(end)/1e6;
                    else
                        tmp=V1; tmp(end)=tmp(end)-(PX(end)-DD(D)*1e6);
                        fill(1+0.4*[-1 1 1 -1], tmp(end+[-1 -1 0 0])/1e6 ,'r','FaceColor',C(t1,:)); hold on; drawnow;
                        Mx1a=tmp(end)/1e6;
                    end
                end
                v1=t1;
            else
                X=[sum(pRPPa(1:v1))+sum(pRPPa(1:v2)) sum(pRPPa(1:v1))+sum(pRPPa(1:t2))]/1e6;
                PX=[PX sum(pRPPa(1:v1))+sum(pRPPa(1:t2))]; V1=[V1 sum(pRPPa(1:v1))]; V2=[V2 sum(pRPPa(1:t2))];
                if PX(end-1)<DD(D)*1e6  % below the limit so plot
                    if PX(end)<DD(D)*1e6  % won't exceded the limit
                        fill(4+0.4*[-1 1 1 -1],V2(end+[-1 -1 0 0])/1e6,'r','FaceColor',C(t2,:)); hold on; drawnow;
                        Mx2a=V2(end)/1e6;
                    else
                        tmp=V2; tmp(end)=tmp(end)-(PX(end)-DD(D)*1e6);
                        fill(4+0.4*[-1 1 1 -1], tmp(end+[-1 -1 0 0])/1e6,'r','FaceColor',C(t2,:)); hold on; drawnow;
                         Mx2a=tmp(end)/1e6;
                    end
                end
                v2=t2;
            end
            axis([0 6 0 22]);
        
        end
    
        % Plot properties
        set(gca,'XTick',[1.5 4.5],'XTickLabel',{'Dose 1','Dose 2'});
        
        text(1,Mx1a,' Age','Rotation',90);
        text(4,Mx2a,' Age','Rotation',90);
        text(2,Mx1,' PG','Rotation',90);
        text(5,Mx2,' PG','Rotation',90);
        title([num2str(DD(D)) 'M']);
    end
    
end

% Plot properties
subplot(3,7,14);
for PG=9:-1:1
    fill([0 1 1 0],PG+0.4*[-1 -1 1 1],'r','FaceColor',C(PG,:)); hold on
    text(1.5,PG,['Priority Group ' num2str(PG)],'VerticalAlignment','middle','FontSize',12);
end
axis([0 5 0 10]);
axis off

print -dpng Bar_Graphs.png
