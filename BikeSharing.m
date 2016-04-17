function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = BikeSharing(x, runlength, seed, ~)
% function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = BikeSharing(x, runlength, seed, other);
% x is a vector containing number of bikes at each station
% runlength is the number of hours of simulated time to simulate
% seed is the index of the substreams to use (integer >= 1)
% other is not used
% Returns
% Average daily cost calculated via expression in writeup
% Standard deviation of daily cost

%   *************************************************************
%   ***             Written by Danielle Lertola               ***
%   ***          dcl96@cornell.edu    June 25th, 2012         ***
%   ***               Written by Bryan Chong                  ***
%   ***        bhc34@cornell.edu    October 29th, 2014        ***
%   *************************************************************

FnGrad = NaN;
FnGradCov = NaN;
constraint = NaN;
ConstraintCov = NaN;
ConstraintGrad = NaN;
ConstraintGradCov = NaN;

if (runlength <= 0) || (round(runlength) ~= runlength) || (seed <= 0) || (round(seed) ~= seed)
    fprintf('runlength and seed should be positive integers\n');
    fn = NaN;
    FnVar = NaN;
else % main simulation
    %% *********************PARAMETERS*********************
    h=16*60;                    %length of simulation (in MINUTES)
    
    %Data will contain
    data=csvread('BikeSharingData.csv',2,1,[2,1,226,229]);
    [m,n]=size(data);
    stations=m;                 %number of stations
    nBikes=3200;                %number of bikes
    nRuns=runlength;            %number of days to be simulated
    capacity=data(:,3);         %capacity of station i
    arrivalRate=data(:,4);      %arrivalRate to station i / hr
    transMat=data(:,5:n);       %P(i,j) of rider biking from station i to j
    start=x.';                  %starting bicycles at each station- column vector
    pEmpty=50/60;               %penalty in $/min for empty rack
    pFull=50/60;                %penalty in $/min for full rack
    r=5;                        %$/bike/km for redistribution
    
    if sum(start)~= nBikes,
        fprintf('starting solution must distribute all bikes');
        fn = NaN;
        FnVar = NaN;
    else
        
        %calculation manhatten distance between two locations
        xlocation=data(:,1);
        ylocation=data(:,2);
        
        xrows=ones(stations,1)*xlocation';
        xcolumns=xlocation*ones(1,stations);
        xdist=abs(xcolumns-xrows);
        
        yrows=ones(stations,1)*ylocation.';
        ycolumns=ylocation*ones(1,stations);
        ydist=abs(ycolumns-yrows);
        
        distance=xdist+ydist;
        
        %calculate Shape and Scale params for time of borrow to return RV's
        rideShape=ones(stations,stations)*(25/20) + (49/45 - 25/20)*eye(stations);
        rideScale=(20^2/25)*distance + (45^2/49)*eye(stations);
        
        %% GENERATE RANDOM NUMBER STREAMS
        % Generate new streams for
        [BikersStream, ArrivalTimeStream, DestinationStream, RideTimeStream] = RandStream.create('mrg32k3a', 'NumStreams', 4);
        % Set the substream to the "seed"
        BikersStream.Substream = seed;
        ArrivalTimeStream.Substream = seed;
        DestinationStream.Substream = seed;
        RideTimeStream.Substream = seed;
        
        %% Generate random # of arrivals data
        OldStream = RandStream.setGlobalStream(BikersStream); % Temporarily store old stream
        % Generate biker turnout per precinct (total arrivals in a day)
        % param for poisson (total arrivals, stations x nRuns)
        lambda = arrivalRate*h/60*ones(1,nRuns);
        Poisson=poissrnd(lambda);
        maxCount= max(max(Poisson));
        ArrivalTime=zeros(stations,maxCount,nRuns);
        
        %% Generate bikers' arrival times.
        RandStream.setGlobalStream(ArrivalTimeStream);
        for i=1:nRuns
            for j=1:stations
                ArrivalTime(j,:,i)=[sort(unifrnd(0, h, 1, Poisson(j,i))) zeros(1,maxCount-Poisson(j,i))];
            end
        end
        
        %% Generate random destination for biker
        RandStream.setGlobalStream(DestinationStream);
        Destination=zeros(size(ArrivalTime));
        for j=1:stations
            cdf=cumsum(transMat(j,:));
            for i=1:nRuns
                x=rand(Poisson(j,i),1);
                for k=1:length(x)
                    sDest=sum(cdf<=x(k))+1;
                    Destination(j,k,i)=sDest;
                end
            end
        end
        %% Generate random ride times (time from borrow to return)
        RandStream.setGlobalStream(RideTimeStream);
        RideTime=zeros(size(Destination));
        
        for j=1:stations
            for i=1:nRuns
                for k=1:Poisson(j,i)
                    RideTime(j,k,i)=gamrnd(rideShape(j,Destination(j,k,i)),rideScale(j,Destination(j,k,i)));
                end
            end
        end
        
        % Restore old random number stream
        RandStream.setGlobalStream(OldStream);
        
        
        %% RUN SIMULATION
        trialCost=zeros(1,nRuns);
        for i=1:nRuns
            % Events columns(type, station, time)
            % type 1:pick up bike from station
            % type 2:return bike to station
            Events=[];
            bikesAvail=start;
            
            clock=0;
            %counter for arrival number of each station
            c=ones(stations,1);
            
            %Indicator counter for empty bike racks (total time on, indicator on/off (1/0), when indicator was last
            %turned on)
            e=zeros(stations,3);
            
            %Indicator counter for full bike racks (total time on, indicator on/off (1/0), when indicator was last
            %turned on)
            f=zeros(stations,3);
            
            %Length of queue of bikes for full bike racks
            q=zeros(stations,1);
            
            for m=1:stations
                if bikesAvail(m)==0 && e(m,2)==0,
                    e(m,2:3)=[1 clock];
                end
                if bikesAvail(m)== capacity(m) && f(m,2)==0,
                    f(m,2:3)=[1 clock];
                end
            end
            %add first arrivals to event list
            for j=1:stations
                if c(j)<=Poisson(j,i),
                    Events=[Events; 1 j ArrivalTime(j,c(j),i)];
                end
            end
            
            Events=sortrows(Events,3);
            while isempty(Events)==0,
                clock=Events(1,3);
                s=Events(1,2);
                type=Events(1,1);
                %event is pick up bike
                if type==1,
                    %there is a bike at station(biker not lost)
                    if bikesAvail(s)~=0,
                        % if the queue for this bike station is empty
                        if q(s) == 0
                            %this rider is picking up last bike, rack is then empty
                            if bikesAvail(s)==1,
                                %start empty rack timer
                                e(s,2:3)=[1 clock];
                            elseif bikesAvail(s) == capacity(s) %rack was full
                                %stop full rack timer, total
                                total=f(s,1) + clock - f(s,3);
                                f(s,:)=[total 0 clock];
                            end
                            bikesAvail(s)=bikesAvail(s)-1;
                        else % there are bikes in the queue, and so rack remains full
                            q(s) = q(s) - 1;
                        end
                        % generate bike drop-off event
                        Events=[Events; 2 Destination(s,c(s),i) clock+RideTime(s,c(s),i)];
                    end
                    % whether or not customer is lost, generate new arrival event
                    c(s)=c(s)+1;
                    if c(s)<=Poisson(s,i),
                        Events=[Events; 1 s ArrivalTime(s,c(s),i)];
                    end
                else %event is bike drop-off
                    
                    %bike rack is already full
                    if bikesAvail(s) == capacity(s)
                        q(s) = q(s) + 1;
                    else
                        %bike rack just became full
                        if bikesAvail(s)+1==capacity(s),
                            %start full rack timer
                            f(s,2:3)=[1, clock];
                        elseif bikesAvail(s)==0, %bike rack was empty
                            %stop empty rack timer, total
                            total=e(s,1) +clock-e(s,3);
                            e(s,:)=[total 0 clock];
                        end
                        bikesAvail(s)=bikesAvail(s)+1;
                    end
                end
                Events(1,:)=[];
                Events=sortrows(Events,3);
            end
            
            %Redistribution
            
            over=max(q+bikesAvail-start,0);
            under=max(start-bikesAvail,0);
            distCost=0;
            for o=1:stations
                %excess bikes at station o
                if over(o)>0,
                    u=1;
                    excess=over(o);
                    while excess>0,
                        if under(u)>0,
                            if excess>=under(u),
                                %There are more extra bikes at station o then
                                %are needed at station u
                                excess=excess-under(u);
                                distCost=distCost + distance(o,u)*under(u);
                                under(u)=0;
                            else
                                %There are more bikes needed at station u
                                %then there are extra bikes at station o
                                under(u)=under(u)-excess;
                                distCost=distCost + distance(o,u)*excess;
                                excess=0;
                            end
                        end
                        u=u+1;
                    end
                end
            end
            %cost matrix for each day
            trialCost(i)=sum(e(:,1))*pEmpty + sum(f(:,1))*pFull + distCost*r;
        end
        
        %Average and Variance of daily cost of bike-sharing program.
        fn=mean(trialCost);
        FnVar=var(trialCost)/nRuns;
        
    end
end
end
