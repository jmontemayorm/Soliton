%% Non-Lineal Schrödinger Equation (NLSE)
%
%   dP   1 d^2P   1 d^2P
%  i-- + - ---- + - ---- + |P|^2 P + V P = 0
%   dz   2 dx^2   2 dy^2
%

%% Settings
% Genetic settings
finalGeneration = 500;
populationSize = 30; % Must be an even number, better if multiple of cores
polynomialTerms = 5;
mutationSize = 10;

% Save and continue settings
saveBestOfGen = 1;
continueProgress = 1;
sinceGen = 0; % Number of save generation, 0 for last found
outF = getOutputFolder(mfilename('fullpath'));

% Parallel settings
maxWorkers = 6;

% Space settings
numOfPoints = 256;
limitZ = 0.3;
windowSize = 10;

%% Calculated settings
x = linspace(-windowSize,windowSize,numOfPoints);
dx = x(2)-x(1);

[xx,yy] = meshgrid(x,x); % Base grid

[xx_0,yy_0] = meshgrid(x-2,x);
r_0 = sqrt(xx_0.^2+yy_0.^2);

[xx_1,yy_1] = meshgrid(x+2,x);
r_1 = sqrt(xx_1.^2+yy_1.^2);

dz = dx^2/8;
doubleSteps = ceil(length(0:dz:limitZ)/2);

% Mode, initial amplitude and omega_0
n = 0;
m = 0;
A0 = 1;
w0 = 1;
[theta_0,~] = cart2pol(xx_0,yy_0);
[theta_1,~] = cart2pol(xx_1,yy_1);

% Psi base
fprintf('Calculating initial PSI state... ');
rw0_0 = (r_0/w0);
rw0_1 = (r_1/w0);
psi_0 = A0 * (rw0_0).^m .* laguerreL(n,2*rw0_0.^2).^m .* exp(-rw0_0.^2 + 1i*(m*theta_0));
psi_1 = A0 * (rw0_1).^m .* laguerreL(n,2*rw0_1.^2).^m .* exp(-rw0_1.^2 + 1i*(m*theta_1));
psi_base =  psi_0 + psi_1;
fprintf('Done!\n');

% Initial potential guess
fprintf('Calculating base potential... ');
V0 = -(abs(psi_base).^2) + 4*n + 2*m + 2 - 2*r_0.^2 - 2*r_1.^2;
fprintf('Done!\n');

%% Initial population
fprintf('Initializing population... ');
% TODO: Get save, if not found, start empty
resetPopulation = false;
if continueProgress == 1
    % look for V_terms, if not found or size differs, resetPopulation true
else
    resetPopulation = true;
end
V_terms = mutationSize * (rand(populationSize,polynomialTerms) - 0.5);

% Initial potentials
V = cell(1,populationSize);
parfor (specimen = 1:populationSize, maxWorkers)
    tempV = V0;
    for p = 1:polynomialTerms
        tempV = tempV + V_terms(specimen,p) * xx.^(p-1);
    end
    V{specimen} = tempV;
end
fprintf('Done!\n');

%% Evolution
fprintf('Starting evolution process...\n\n');
tic
psi_results = cell(1,populationSize);
specimen_fitness = zeros(1,populationSize);
populationIndices = 1:populationSize;
for generation = 1:finalGeneration
    % Parallel propagation
    fprintf('Propagating generation number %04i... ',generation);
    parfor (specimen = 1:populationSize, maxWorkers)
        psi_results{specimen} = SolitonPropagator2D(psi_base,V{specimen},windowSize,dz,doubleSteps)
        specimen_fitness(specimen) = sum(sum(abs(abs(psi_results{specimen}) - abs(psi_base))));
    end
    fprintf('Finished soliton propagation.\n');
    
    % Evaluate
    [~,bestIdx] = sort(specimen_fitness);
    fprintf('\tBest fit in this generation is: %0.7f\n',specimen_fitness(bestIdx(1)));
    
    % Save the best (in file)
    if saveBestOfGen == 1
        bestOfGen = V_terms(bestIdx(1),:);
        save(fullfile(outF,sprintf('BestOfGen_%05i.mat',generation)),'bestOfGen');
    end
    
    % Survival of the fittest
    fprintf('\tKilling the unfit... ');
    killed = 0;
    killIdx = false(1,populationSize);
    while killed < populationSize / 2
        % Go from unfittest to fittest
        for specimen = flip(bestIdx)
            % Always include probability of survival (also for unfittest)
            if  ~killIdx(specimen) && exp(find(specimen == bestIdx,1)/populationSize - 1.1) >= rand
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
    
    % Kill and substitute via reproduction
    fprintf('Replacing corpses with new offspring... ');
    replaceWithBaby = find(killIdx);
    for newBaby = 1:length(replaceWithBaby)
        % Only search in the ones not to be killed
        allowedIdx = ~killIdx;
        
        % Get first parent
        firstParentIdx = 0;
        lookingForFirstParent = true;
        while lookingForFirstParent
            for candidate = populationIndices(allowedIdx)
                if rand < 0.4 % Equal probability, give a bias to fittest via order?
                    lookingForFirstParent = false;
                    firstParentIdx = candidate;
                    allowedIdx(firstParentIdx) = false;
                    break
                end
            end
        end
        
        % Get second parent
        secondParentIdx = 0;
        lookingForSecondParent = true;
        while lookingForSecondParent
            for candidate = populationIndices(allowedIdx)
                if rand < 0.4 % Equal probability, give a bias to fittest via order?
                    lookingForSecondParent = false;
                    secondParentIdx = candidate;
                    break
                end
            end
        end
        
        % Make baby
        firstParent = randi([0, 1],[1, polynomialTerms]);
        secondParent = ~firstParent;
        V_terms(replaceWithBaby(newBaby),:) = firstParent .* V_terms(firstParentIdx,:) + secondParent .* V_terms(secondParentIdx,:);
    end
    fprintf('Done!\n');
    
    % Generate and apply mutations
    fprintf('\tApplying mutations and generating new potentials... ');
    V_terms_mutations = mutationSize * (rand(populationSize,polynomialTerms) - 0.5);
    mutate = rand(populationSize,polynomialTerms) < 0.1;
    V_terms = V_terms + mutate .* V_terms_mutations;
    % esquema de enfriamiento
    % Suma / resta de valores rand
    
    % Calculate new potentials
    parfor (specimen = 1:populationSize, maxWorkers)
        tempV = V0;
        for p = 1:polynomialTerms
            tempV = tempV + V_terms(specimen,p) * xx.^(p-1);
        end
        V{specimen} = tempV;
    end
    fprintf('Done!\n\n');
end
toc

fprintf('Evolution sequence complete. Achieved generation %04i!\n',generation);