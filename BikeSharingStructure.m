function [minmax d m VarNature VarBds FnGradAvail NumConstraintGradAvail, StartingSol, budget, ObjBd, OptimalSol] = BikeSharingStructure(NumStartingSol, seed)
% Inputs:
%	a) NumStartingSol: Number of starting solutions required. Integer, >= 0
%	b) seed: Seed for generating random starting solutions. Integer, >= 1
% Return structural information on optimization problem
%     a) minmax: -1 to minimize objective , +1 to maximize objective
%     b) d: positive integer giving the dimension d of the domain
%     c) m: nonnegative integer giving the number of constraints. All
%        constraints must be inequality constraints of the form LHS >= 0.
%        If problem is unconstrained (beyond variable bounds) then should be 0.
%     d) VarNature: a d-dimensional column vector indicating the nature of
%        each variable - real (0), integer (1), or categorical (2).
%     e) VarBds: A d-by-2 matrix, the ith row of which gives lower and
%        upper bounds on the ith variable, which can be -inf, +inf or any
%        real number for real or integer variables. Categorical variables
%        are assumed to take integer values including the lower and upper
%        bound endpoints. Thus, for 3 categories, the lower and upper
%        bounds could be 1,3 or 2, 4, etc.
%     f) FnGradAvail: Equals 1 if gradient of function values are
%        available, and 0 otherwise.
%     g) NumConstraintGradAvail: Gives the number of constraints for which
%        gradients of the LHS values are available. If positive, then those
%        constraints come first in the vector of constraints.
%     h) StartingSol: One starting solution in each row, or NaN if NumStartingSol=0.
%        Solutions generated as per problem writeup
%     i) budget: Column vector of suggested budgets, or NaN if none suggested
%     j) ObjBd is a bound (upper bound for maximization problems, lower
%        bound for minimization problems) on the optimal solution value, or
%        NaN if no such bound is known.
%     k) OptimalSol is a d dimensional column vector giving an optimal
%        solution if known, and it equals NaN if no optimal solution is known.

%   *************************************************************
%   ***             Written by Danielle Lertola               ***
%   ***          dcl96@cornell.edu    June 25th, 2012         ***
%   ***               Written by Bryan Chong                  ***
%   ***       bhc34@cornell.edu    October 29th, 2014         ***
%   *************************************************************

% Number of stations is user selectable, but in this version
% we use 225.

data=csvread('BikeSharingData.csv',2,1,[2,1,226,229]);
[m,~]=size(data);
stations=m;


minmax = -1; % minimize total cost of penalties for empty/full racks and reallocation
d = stations; % # bikes allocated to station i in i=[1,stations]
m = 3; % final station allocation must fall in [0, c(225)], and variables must sum to total # of bikes
VarNature = ones(d, 1); % integer variables
VarBds = [zeros(d,1) data(1:d,3)]; % capacity of stations 
FnGradAvail = 0; % No derivatives 
NumConstraintGradAvail = 0; % No constraint gradients
%budget = [100; 1000; 5000];
budget = [5000];
ObjBd=NaN;
OptimalSol=NaN;

if (NumStartingSol < 0) || (NumStartingSol ~= round(NumStartingSol)) || (seed <= 0) || (round(seed) ~= seed),
    fprintf('NumStartingSol should be integer >= 0, seed must be a positive integer\n');
    StartingSol = NaN;
    
else
    bikes=3200; %number of bikes
    cap=data(:,3); %station capacity
       
    if (NumStartingSol == 0),
        StartingSol = NaN;
    else
        %Initial allocation z(i) according to write up
        StartingSol=ones(NumStartingSol,1)*min(floor(bikes/stations),cap.');
        if(bikes-sum(StartingSol(1,:))> 0),
            OurStream = RandStream.create('mlfg6331_64'); 
            % Use a different generator from simulation code to avoid stream clashes
            OurStream.Substream = seed;
            OldStream = RandStream.setGlobalStream(OurStream);      
            
            %Distribute remaining bikes
            remain=bikes-sum(StartingSol(1,:));
            %for i=1:NumStartingSol
            %    x=rand(remain,1);
            %    for j=1:remain,
            %        y=cap.'-StartingSol(i,:);
            %        yProb=y./sum(y);
            %        cdf=cumsum(yProb);
            %        sAddTo=sum(cdf<=x(j))+1;
            %        StartingSol(i,sAddTo)=StartingSol(i,sAddTo)+1;
            %    end
            
            space = 0;
            for i = 1:stations
                space = space + cap(i,:)-StartingSol(:,i);
            end
            for i = 1: remain
                loc = unidrnd(space);
                a = 0;
                numrack = 1;
                while (a <= loc) & (numrack < 266)
                    topSpaceRack = a + cap(numrack) - StartingSol(numrack);
                    a = topSpaceRack;
                    numrack = numrack + 1;       
                end
                StartingSol(numrack-1) = StartingSol(numrack-1) + 1;
                space = space -1;
            end
            RandStream.setGlobalStream(OldStream); % Restore previous stream
        end  
        
    end %if NumStartingSol 
    
end %if inputs ok
