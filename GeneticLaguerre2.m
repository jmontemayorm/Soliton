%% Non-Lineal Schrödinger Equation (NLSE)
%
%   dP   1 d^2P   1 d^2P
%  i-- + - ---- + - ---- + |P|^2 P + V P = 0
%   dz   2 dx^2   2 dy^2
%

%% Settings
% Soliton parameters' limits
A_lowerBound = 0.1;
A_upperBound = 100;
A_mutation = 0.01;
B_lowerBound = 0.1;
B_upperBound = 20;
B_mutation = 0.01;

% Function settings
% TODO: Refactor
lambda = 0.5;
V = 0;
N = 1;
s = 0.02;
D = 1;

% Space settings
windowSize = 20;
numOfPoints = 256;

% Population
populationSize = 100;

% Mutation
paternalProbability = 0.6;
mutationProbability = 0.05;

% Generations
enableMaxGenerations = 1;
maxGenerations = 2000;

%% Calculated settings
% Space
x = linspace(-windowSize/2,windowSize/2,numOfPoints);
dx = x(2) - x(1);

[xx,yy] = meshgrid(x,x);
rr = xx.^2 + yy.^2;

% TODO: inline function
% Function handle (anonymous)
f = @(u) -(lambda * u) + (D * del2(u,dx)) + (N * (u .* abs(u).^2) ./ (1 + s * abs(u).^2)) + (V * u);

%% Initial population
theLiving = cell(populationSize,2);
for specimen = 1:populationSize
    theLiving{specimen}(1) = A_lowerBound + (A_upperBound - A_lowerBound) * rand;
    theLiving{specimen}(2) = B_lowerBound + (B_upperBound - B_lowerBound) * rand;
end

specimenFitness = zeros(populationSize,1);

%% Evolution
tic

generation = 0;
bestFitness = 0;

while true
    generation = generation + 1;
    
    % % % Evaluation % % %
    for specimen = 1:populationSize
        % One-time evaluation and storage is faster
        fU = f(theLiving{specimen}(1) .* exp(-rr ./ theLiving{specimen}(2).^2));
        specimenFitness(specimen) = 1 ./ sum(sum( conj(fU) .* fU ));
    end
    
    % Sort fitness
    [~,bestIdx] = sort(specimenFitness,'descend');
    
    % TODO
    % Check and save (RAM) the all-time best
    if specimenFitness(bestIdx(1)) > bestFitness
        bestFitness = specimenFitness(bestIdx(1))
        bestSpecimen = theLiving{bestIdx(1)};
        bestGeneration = generation;
    end
    
    % % % Survival of the fittest % % %
    % Acquire targets
    killed = 0;
    killIdx = false(1,populationSize);
    while killed < populationSize / 2
        % Go from unfittest to fittest
        unfitToFit = flip(bestIdx); % TODO: Only the non elite
        
        for s = 1:(populationSize)% TODO: - eliteAmount)
            specimen = unfitToFit(s);
            % Always include probability of survival (also for unfittest, elitism exception)
            if (killIdx(specimen) == false) && (exp(find(specimen == bestIdx,1)/populationSize - 1.1) >= rand)
                % Acquire target
                killIdx(specimen) = true;
                
                % Increase counter and check
                killed = killed + 1;
                if killed >= populationSize / 2
                    break
                end
            end
        end
    end
    
    % Make clones of the best
    % TODO: Make this flexible
    % Make 2 clones of the best
    cloneIdx = find(killIdx,2);
    theLiving{cloneIdx(1)} = bestSpecimen;
    theLiving{cloneIdx(2)} = bestSpecimen;
    killIdx(cloneIdx(:)) = false;
    
    % Kill and substitute via reproduction
    replaceWithBaby = find(killIdx);
    for newBaby = 1:length(replaceWithBaby)
        % Only search in the ones not to be killed
        allowedIdx = ~killIdx;
        
        % Get first parent
        firstParentIdx = 0;
        lookingForFirstParent = true;
        while lookingForFirstParent
            % The most fit are the firsts in line
            orderedCandidatesIdx = bestIdx(allowedIdx);
            
            for candidate = 1:length(orderedCandidatesIdx)                
                if rand < paternalProbability
                    lookingForFirstParent = false;
                    firstParentIdx = orderedCandidatesIdx(candidate);
                    allowedIdx(firstParentIdx) = false;
                    break
                end
            end
        end
        
        % Get second parent
        secondParentIdx = 0;
        lookingForSecondParent = true;
        while lookingForSecondParent
            % The most fit are the firsts in line
            orderedCandidatesIdx = bestIdx(allowedIdx);
            
            for candidate = 1:length(orderedCandidatesIdx)
                if rand < paternalProbability
                    lookingForSecondParent = false;
                    secondParentIdx = orderedCandidatesIdx(candidate);
                    break
                end
            end
        end
        
        % Make baby
        firstParent = randi([0, 1], [1, 2]);
        secondParent = ~firstParent;
        theLiving{replaceWithBaby(newBaby)} = firstParent .* theLiving{firstParentIdx} + secondParent .* theLiving{secondParentIdx};
    end
    
    % % % Mutations % % %
    % Mutate
    for specimen = 1:populationSize
        mutate = rand(1, 2) < mutationProbability;
        if mutate(1)
            ss = 2 * randi([0,1],[1,1]) - 1;
            theLiving{specimen}(1) = theLiving{specimen}(1) + ss * A_mutation;
            if theLiving{specimen}(1) < A_lowerBound
                theLiving{specimen}(1) = A_lowerBound;
            elseif theLiving{specimen}(1) > A_upperBound
                theLiving{specimen}(1) = A_upperBound;
            end
        end
        
        if mutate(2)
            ss = 2 * randi([0,1],[1,1]) - 1;
            theLiving{specimen}(2) = theLiving{specimen}(2) + ss * B_mutation;
            if theLiving{specimen}(2) < B_lowerBound
                theLiving{specimen}(2) = B_lowerBound;
            elseif theLiving{specimen}(2) > B_upperBound
                theLiving{specimen}(2) = B_upperBound;
            end
        end
    end
    
    % % % Breaking mechanisms % % %
    if enableMaxGenerations == 1 && generation == maxGenerations
        break
    end
end
toc