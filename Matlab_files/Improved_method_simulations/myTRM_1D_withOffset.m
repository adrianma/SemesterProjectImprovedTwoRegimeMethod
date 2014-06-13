function [rMytest2,vcenters,vVoxelN,vMeanPosition]= ...
    myTRM_1D_withOffset(vX,vActiveP,vBoundaries,vTimes,OffsetDecayTuning)
%
%    Improved two-regime method (1D)
%
%    vx  : Position in space of all particles (x position).
%    vActiveP: Set to boolean(1) if particle is in OmegaM (active).
%    vBoundaries : vector with the boundaries of the 1D problem.
%    vTimes : vector with the times defining the simulation.
%    OffsetDecayTuning : parameter for the offset decay rate constant. 
%
%    rMytest2 : data struct containing useful stored data for mytest2.
%    vcenters : voxel centers.
%    vVoxelN : number of particles in every voxel.
%    vMeanPosition : mean position for offset analysis.
%

%% Constants and definitions
% Boundary between the "forbidden region" and OmegaC
sxBoundaryFR1 = vBoundaries(1);
% Boundary between OmegaC and OmegaM
sxBoundaryCM = vBoundaries(2);
% Boundary between OmegaM and the "forbidden region"
sxBoundaryFR2 = vBoundaries(3);
% Diffusion coefficient D = 1 [length_unit^2/time_unit]
sD = 1.0;
% Duration of the simulation
st = vTimes(1); stEnd = vTimes(2); st0 = vTimes(3); sDeltaT = vTimes(4); 
vT = vTimes(5);
% Gillespie length of voxel side and minimum virtual voxel length
sh = 1.0; shMin = sh*0.01;
% Number of particles in the regions OmegaM and OmegaC
sNParticles = length(vX);
% Number of compartments (in OmegaC AND OmegaM).
sNK = (sxBoundaryFR2-sxBoundaryFR1)/sh;

%% Compute useful vectors and matrices for properties about particles
% vIcisGillespie: boolean indicating if the voxel is in the Gillespie
%                 domain.
%           1: the particle is in OmegaC (is Gillespie)
%           0: the particle is in OmegaM (is Brownian)
vIcisGillespie = boolean(zeros(1,sNK));
vIx = find((sxBoundaryFR1:sh:(sxBoundaryFR2-1))<sxBoundaryCM);
vIcisGillespie(vIx) = boolean(1);

% mN: matrix showing how many particles of every type in every
%     voxel of the Gillespie-mesh.
mN = zeros(sNK,1);
for i = 1:sNParticles
    if(vIcisGillespie(floor(vX(i))))
        mN(floor(vX(i)),1) = mN(floor(vX(i)),1) + 1;
        
        % If particle in Gillespie domain, set to inactive.
        vActiveP(i) = boolean(0);
    end
end

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

% vNConnects: Total number of possible connections for a voxel.
vNConnects = mIcConnects(:,1) + mIcConnects(:,2);

%% Voxel offset 'mOffset' and adaptive mesh size 'sh_'
% Adapted array size to hold t0 as well.
vOffset = zeros(sNK,2); 
sh_ = sh.*ones(size(vOffset,1),size(vOffset,2)) - 2*vOffset(1);

% Offset decay rate constant.
% (adjusted to relaxation simulation tests.)
sOffset_k = sD/sh^2 * OffsetDecayTuning; 

%% Calculate putative times
% Assume regular lattice, the matrix containing the propensity for a 
% diffusive event is also computed.
skJump = sD/sh^2;
mAlphaEpsilon = skJump.*mN;

% Putative times matrix. 
% Entry(1): next time.
% Entry(2:3): weight for directions
vtEpsilon = zeros(sNK,1);
vtEpsilonW = zeros(sNK,2); 

for j = 1:sNK
    vtEpsilon(j,1) = st + ...
        1.0/(mAlphaEpsilon(j,1)*vNConnects(j))*exprnd(1,1);
    vtEpsilonW(j,1) = mIcConnects(j,1);
    vtEpsilonW(j,2) = mIcConnects(j,2);
end

% stM : time for the next M-event (Particle).
% stC : time for the next C-event (Gillespie).
stM = sDeltaT;
stC = min(vtEpsilon);

rData = struct('vX',vX,'vActiveP',vActiveP,'sNK',sNK,'sNParticles',...
    sNParticles);
rGillespie = struct('vIcisGillespie',vIcisGillespie,'mIcConnects',...
    mIcConnects,'vNConnects',vNConnects,'mN',mN,'skJump',skJump,...
    'mAlphaEpsilon',mAlphaEpsilon,'vtEpsilon',vtEpsilon,...
    'vtEpsilonW',vtEpsilonW,...
    'sxBoundaryFR1',sxBoundaryFR1,'sxBoundaryCM',sxBoundaryCM,...
    'sxBoundaryFR2',sxBoundaryFR2,'sDeltaT',sDeltaT,'sh',sh,...
    'sD',sD,'shMin',shMin,'st',st,'stEnd',stEnd,'stM',stM,...
    'stC',stC,'sAdjacent',sAdjacent,'st0',st0,'sh_',sh_,'vOffset',...
    vOffset,'sOffset_k',sOffset_k);

%% Simulation
ix = 0;
bins = vBoundaries(1):rGillespie.sh:vBoundaries(3);

[px_particle,X] = hist(rData.vX,bins);

while(rGillespie.st<rGillespie.stEnd)
    
    % Next C-event occurs.
    if(rGillespie.stC <= rGillespie.stM)
        
        % Iterations.
        while(rGillespie.stC<rGillespie.stM)
                       
            % Update current time
            rGillespie.st = rGillespie.stC;
            
            [rDataN,rGillespieN] = ...
                GillespieAlgorithm_1D(rData,rGillespie);
            
            % Recompute time for the next C-event and save data.
            rGillespieN.stC = min(rGillespieN.vtEpsilon);
            
            rData = rDataN;
            rGillespie = rGillespieN;
        end
        
    % Next M-event occurs.
    else 
        
        % Update current time
        rGillespie.st = rGillespie.stM;
        fprintf('Particle Diffusion to t=%f \n',rGillespie.st);
        
        [rDataN,rGillespieN] =...
            ParticleAlgorithm_1D(rData,rGillespie);
        
        rGillespieN.stM = rGillespieN.stM + rGillespieN.sDeltaT;
        
        rGillespieN.stC = min(rGillespieN.vtEpsilon);
        
        rData = rDataN;
        rGillespie = rGillespieN;
        
        ix = ix + 1;
        
    end
    
    ix = max(ix,1);
    
    % Compute the mean_position(t) offset (which is vMeanPosition(:,1)) for
    % comparison purposes.
    [px_particle] = histc(rData.vX(rData.vActiveP),bins);
    px_particle1 = px_particle(1:size(px_particle,1)-1);    
    
    vVoxelN = px_particle1 + rGillespie.mN;
    vcenters = bins + 0.5;
    vcenters = vcenters(1:length(vcenters)-1);
    
    d = sum(vVoxelN);
    e = vcenters.*vVoxelN';
    e = sum(e);
    vMeanPosition(ix,1) = e/d;
    vMeanPosition(ix,2) = rGillespie.st;
    
end

% Store values in the struct for analysis in mytest2.
vXNew = rData.vX;
vActivePNew = rData.vActiveP;
mmN = rGillespie.mN;
vvIcisGillespie = rGillespie.vIcisGillespie;

rMytest2 = struct('vXNew',vXNew,'vActivePNew',vActivePNew,'mmN',mmN,...
    'vvIcisGillespie',vvIcisGillespie,'vOffset',vOffset);

end