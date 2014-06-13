function [mmN,sParticleRuntime] = ...
    gillespieTest_1D(sPlots,sSimulation,sDistribution);
%
%   Set sDistribution to 0 for random uniform distribution or
%   to 1 for delta distribution.
%
%   Problem is 1D.
%
%   Use [mmN] = gillespieTest_1D(1,0,0); for rand. uni. distr.
%   or  [mmN] = gillespieTest_1D(1,0,1); for delta distr.
%

%% Constants and definitions
% Number of particles.
sNParticles = 10000;
% Starting position (delta distribution).
x0 = 5.0;
% Boundaries to the forbidden regions.
sxBoundaryFR1 = 1.0; sxBoundaryFR2 = 10.0;
% Diffusion coefficient D = 1 [length_unit^2/time_unit]
sD = 1.0;
% Duration of the simulation / Time step / Time vector.
st = 0; stEnd = 30.0; sDeltaT = 0.01; vT = st:sDeltaT:(stEnd-st);
% Steps
sN = (stEnd-st)/sDeltaT;
% Gillespie length of voxel side and minimum virtual voxel length
sh = 1.0; skJump = sD/sh^2;
% Number of compartments.
sNK = (sxBoundaryFR2-sxBoundaryFR1)/sh;
% Initialize output matrix.
mmN = zeros(sNK,sN);
% Number of trials the experiment is run.
sNTrials = 1;
% Bin width.
sbinWidth = 1.0;

% mConnects: connection values to neighboring voxels.
%            The 1. column represents connections to the left.
%            The 2. column represents connections to the right.
%            0 : particle CAN NOT move in this direction.
%            1 : particle CAN move in this direction.
% If particle located in the first voxel of the mesh, set left connection,
%  to 0 to indicate that it only can move to the right (or viceversa).
%       /// m = 2 for 1D; m = 4 for 2D; m = 6 for 3D ///
sAdjacent = 2;
mIcConnects = boolean(ones(sNK,sAdjacent));
mIcConnects(1,1) = boolean(0);
mIcConnects(end,2) = boolean(0);

vIcisGillespie = boolean(ones(sNK,1));

% vNConnects: Total number of possible connections for a voxel.
vNConnects = mIcConnects(:,1) + mIcConnects(:,2);

%% Simulation parameter definitions
bins = sxBoundaryFR1:sh:sxBoundaryFR2;
mN = zeros(sNK,1);

% Delta distribution for the particles
if(sDistribution==1)
    vX = ones(sNParticles,1).*x0;
 
% Uniform distribution for the particles
elseif(sDistribution==0)
    vX = 1+rand(sNParticles,1).*(sxBoundaryFR2-1);
end

for j = 1:sNParticles
    mN(floor(vX(j)),1) = mN(floor(vX(j)),1) + 1;
end

mAlphaEpsilon = skJump.*mN;

% Putative times matrix.
vtEpsilon = zeros(sNK,1);
for j = 1:sNK
    vtEpsilon(j,1) = st + ...
        1.0/(mAlphaEpsilon(j,1)*mIcConnects(j,1))*exprnd(1,1);
end

%skJump = sqrt(2*sD*rGillespie.st);
skJump = sD/sh^2;

rData = struct('sNParticles',sNParticles);
rGillespie = struct('sxBoundaryFR1',sxBoundaryFR1,'sxBoundaryFR2',...
    sxBoundaryFR2,'sDeltaT',sDeltaT,'sD',sD,'st',st,'stEnd',stEnd,...
    'sh',sh,'mIcConnects',mIcConnects,'vNConnects',vNConnects,...
    'mN',mN,'mAlphaEpsilon',mAlphaEpsilon,'vtEpsilon',vtEpsilon,...
    'skJump',skJump,'vIcisGillespie',vIcisGillespie);

%% Simulation
ix = 1;
tic;

% Plot the particle concentrations as histogram plot.
if(sPlots==1 && sSimulation==1)
    hist(rGillespie.mN,bins); grid on;
    xlabel('x'); ylabel('Number of particles');
    title('Experiment ');
    
    if(sSimulation==1)
        pause(0.3);
    end
end

while(rGillespie.st<rGillespie.stEnd)
    
    mmN(:,ix) = rGillespie.mN;
    
    [rDataN,rGillespieN]=GillespieAlgorithm_1D(rData,rGillespie);
    rGillespieN.st = rGillespieN.st + rGillespieN.sDeltaT;
    
    rData = rDataN;
    rGillespie = rGillespieN;
    
    % Plot the particle concentrations as histogram plot.
    if(sPlots==1 && sSimulation==1)
        figure(1);
        hist(rGillespie.mN,bins); grid on;
        xlabel('x'); ylabel('Number of particles');
        title('Experiment');
        
        pause(0.3);
    end

    ix = ix + 1;
end

mmN(:,ix) = rGillespie.mN;

% Compute the simulation time and print it out.
sParticleRuntime = toc;
if(sSimulation==0)
    fprintf('\n \n');
    fprintf('Particle runtime %f s \n\n',sum(sParticleRuntime));
end

%% Plot histogram
if(sPlots==1 && sSimulation==0)
    figure(3);
    stem(sxBoundaryFR1:sh:(sxBoundaryFR2-1),mmN(:,end),'.','LineWidth',3); 
    grid on;
    
    hXlabel = xlabel('Voxel'); hYlabel = ylabel('Number of particles');
    hLegend = legend([num2str(sum(mN)),' particles']);
    hTitle = title('Simulation');
    
    set([hXlabel,hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
end

end

function [rDataN,rGillespieN]=GillespieAlgorithm_1D(rData,rGillespie)
%
%    Refined Gillespie Algorithm (1D)
%    Adrian Martinez Gomez,  16/September/2013
%    Revised by Michael Klann, 25/September/2013
%

%% Find voxel index for particle with minimum value.

[~,index]=min(reshape(rGillespie.vtEpsilon,numel(rGillespie.vtEpsilon),1));
[sVoxel,~] = ind2sub(size(rGillespie.vtEpsilon),index);

fprintf('... from voxel %i',sVoxel);

%% Particle jumps out from the activated voxel.
[rGillespieNew] =...
    G_Jump_out(sVoxel,rGillespie);
rGillespie = rGillespieNew;

%% Determine Jump direction sDir.
% sDir is:
% 2 : if particle comes from the left (has jumped to the right).
% 1 : if particle comes from the right (has jumped to the left).

if(rGillespie.vNConnects(sVoxel)==1)
    % Only jump to the right/east possible.
    if (rGillespie.mIcConnects(sVoxel,:)==boolean([0,1]) )
        sDir = 2;
    % Only jump to the left/west possible.
    else
        sDir = 1;
    end
else
    % Using precomputed direction weights
    vDirections = [0.5,0.5];
    sDir = selectDirection(vDirections);
end

%% Create new particle / regard the domains

if(sDir==2)
    sVoxelNew = sVoxel+1;
else
    sVoxelNew = sVoxel-1;
end

% Jump to OmegaM
if (~rGillespie.vIcisGillespie(sVoxelNew))
    [rDataNew] = P_Create_particle(sVoxel,sDir,rData);
    rData = rDataNew;
    % Jump to OmegaC
else
    [rGillespieNew] = ...
        G_Jump_in(sVoxelNew,rGillespie);
    rGillespie = rGillespieNew;
end

fprintf('... to voxel %i \n',sVoxelNew);


% Save the (updated) datastructures.
rDataN = rData;
rGillespieN = rGillespie;

end

function [rGillespieNew] = ...
    G_Jump_out(sVoxel,rGillespie)
%
%   /// Description of input & output of the function /// 
%

%% Delete the particle that just abandoned the voxel.
rGillespieNew = rGillespie;

% (3.d.i.1)
rGillespieNew.mN(sVoxel,1) = ...
    rGillespieNew.mN(sVoxel,1) - 1;

%% Redraw the putative (jump) times.
rGillespieNew.vtEpsilon(sVoxel,1) = rGillespieNew.st+1.0/...
    (rGillespieNew.mN(sVoxel,1)*(rGillespieNew.skJump)).*exprnd(1,1);

end

function [rGillespieNew] =...
    G_Jump_in(sVoxel,rGillespie)
%
%   /// Description of input & output of the function /// 
%

%% Add the particle that just entered the voxel.
rGillespieNew = rGillespie;

%(3.c.iv) Update Numbers
rGillespieNew.mN(sVoxel,1) = ...
    rGillespieNew.mN(sVoxel,1) + 1;

%% Redraw the putative (jump) times.
rGillespieNew.vtEpsilon(sVoxel,1) = rGillespieNew.st+1.0/...
    (rGillespieNew.mN(sVoxel,1)*(rGillespieNew.skJump)).*exprnd(1,1);

end

function [rDataNew] = P_Create_particle(sVoxel,sDir,rData)
%
%   /// Description of input & output of the function /// 
%

%% Position of the transfered particle in the OmegaM domain and set active.
rDataNew = rData;

[~,sIx] = min(rDataNew.vActiveP);
rDataNew.vActiveP(sIx) = boolean(1);

if(sDir==2)
    rDataNew.vX(sIx) = sVoxel + 1;
elseif(sDir==1)
    rDataNew.vX(sIx) = sVoxel - 1;
end


end

function [sDir] = selectDirection(vDirections)
%
% Credits: @ Michael Klann
%
% vDirections is an array holding the weight of the different directions;
% vDirections does not have to be normalized

% compute the sum and weight with a uniform random number
% (Gillespie 1976)
random_weight = sum(vDirections)*rand();
% iterate through the directions weights
% until the growing sum maches or exceeds the random_weight
% in order to select each direction with a probability
% proportional to its weight
sum_d = 0;
sDir = length(vDirections); %(preset the direction to the last direction)
for i = 1:length(vDirections)
    sum_d = sum_d + vDirections(i);
    if (sum_d >= random_weight)
        sDir = i;
        return;
    end
end

% if you reach this point, for whatever reason (you should not ...
%   - because rand() <= 1 ==> random_weight <= sum(directions)
%   - ==> final sum_d >= random weight (except for numerical reasons)
% then don't be afraid, sDir is preset to the last direction in order to
% capture that numerical error.
end