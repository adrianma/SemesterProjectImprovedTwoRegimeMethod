function [rDataN,rGillespieN]=GillespieAlgorithm_1D(rData,rGillespie)
%
%    Refined Gillespie Algorithm (1D)
%

%% Find voxel index for particle with minimum value.

[~,index]=min(reshape(rGillespie.vtEpsilon,numel(rGillespie.vtEpsilon),1));
[sVoxel,~] = ind2sub(size(rGillespie.vtEpsilon),index);

%% Particle jumps out from the activated voxel.
[rGillespieNew,sParticleOffsetNew] =...
    G_Jump_out(sVoxel,rGillespie);
rGillespie = rGillespieNew;
sParticleOffset = sParticleOffsetNew;

%% Determine Jump direction sDir.
% sDir is:
% 2 : if particle comes from the left (has jumped to the right).
% 1 : if particle comes from the right (has jumped to the left).

if(rGillespie.vNConnects(sVoxel)==1)
    if (rGillespie.mIcConnects(sVoxel,:)==boolean([0,1]) )
        % Only jump to the right/east possible.
        sDir = 2;
    else
        % Only jump to the left/west possible.
        sDir = 1;
    end
else 
    % Jump to both directions is possible.
    vDirections = [rGillespieNew.vtEpsilonW(sVoxel,1),rGillespieNew.vtEpsilonW(sVoxel,2)]; % MK_correction: using precomputed direction weights
    sDir = selectDirection(vDirections);
end

%% Mirror and cut off the particle offset:

% Mirroring the offset in jump direction.
sParticleOffset = -sParticleOffset; 
if (sDir==2) % particle will jump to the right.
    % if offset > 0, set to 0.
    sParticleOffset = min(sParticleOffset,0.0);
    
    % Particle jumps to the right
    sVoxelNew = sVoxel+1;
    
else % particle will jump to the left.
    % if offset < 0, set to 0.
    sParticleOffset = max(sParticleOffset,0.0);
    
    % Particle jumps to the left
    sVoxelNew = sVoxel-1;
end
% ==> Only keep the offset if you jump into the short direction

%% Create new particle / regard the domains

    if (~rGillespie.vIcisGillespie(sVoxelNew)) % jump to to OmegaM
        [rDataNew] = P_Create_particle(sVoxel,sDir,sParticleOffset,rData);
        rData = rDataNew;
    else % jump to to OmegaC
        [rGillespieNew] = ...
            G_Jump_in(sVoxelNew,rGillespie,sParticleOffset); 
            % MK_Correction: mirror the offset. In 1D it is always
        rGillespie = rGillespieNew;
    end
    
% Save the (updated) datastructures.
rDataN = rData;
rGillespieN = rGillespie;

end

% [x] OK: MK_correction: partially rewritten
function [rGillespieNew,sParticleOffsetNew] = ...
    G_Jump_out(sVoxel,rGillespie)
%
%   /// Description of input & output of the function /// 
%

%% Delete the particle that just abandoned the voxel.
rGillespieNew = rGillespie;

% (3.d.i.1)
rGillespieNew.mN(sVoxel,1) = ...
    rGillespieNew.mN(sVoxel,1) - 1;

if (rGillespieNew.sOffset_k > 0)

% (3.d.i.2) Relax offset to current time.
%   Current time = rGillespie.st
%   time of entry = rGillespieNew.vOffset(sVoxel,2)
%   Offset decay rate constant = rGillespieNew.sOffset_k
rGillespieNew.vOffset(sVoxel,1) = rGillespieNew.vOffset(sVoxel,1) ...
    *exp( -rGillespieNew.sOffset_k ...
* (rGillespieNew.st-rGillespieNew.vOffset(sVoxel,2) ));

% Update time
rGillespieNew.vOffset(sVoxel,2) = rGillespie.st; 
sParticleOffsetNew = rGillespieNew.vOffset(sVoxel,1);

else
    rGillespieNew.vOffset(sVoxel,1) = 0.0;
    sParticleOffsetNew = 0.0;
end

%% Recomputing The Jump Times:
%   Calculating the effective sh_.
%    Positive Offset means smaller sh_, means jump to the right is more
%    likely.
rGillespieNew.sh_(sVoxel) = rGillespieNew.sh - ...
    2*rGillespieNew.vOffset(sVoxel,1);

% Cut-off such that h_>h_min>0.
if(rGillespieNew.sh_(sVoxel)<rGillespieNew.shMin)
    rGillespieNew.sh_(sVoxel) = rGillespieNew.shMin;
else
    if (rGillespieNew.sh_(sVoxel)>2-rGillespieNew.shMin)
        rGillespieNew.sh_(sVoxel) =2-rGillespieNew.shMin;
    end
end

% Local jump rate constant (if connected:
% to the left: 
skJump_1 = rGillespieNew.sD*rGillespieNew.mIcConnects(sVoxel,1)/...
    (2-rGillespieNew.sh_(sVoxel))^2;
% to the right:
skJump_2 = rGillespieNew.sD*rGillespieNew.mIcConnects(sVoxel,2)/...
    (rGillespieNew.sh_(sVoxel))^2;
% As noted above, due to division by sh_, for sh_ < 0.5 jump to the right
% is more likely.

%% Redraw the putative (jump) times.
rGillespieNew.vtEpsilon(sVoxel,1) = rGillespieNew.st+1.0/...
    (rGillespieNew.mN(sVoxel,1)*(skJump_1+skJump_2)).*exprnd(1,1);

% And store the jump weights:
rGillespieNew.vtEpsilonW(sVoxel,1) = skJump_1;
rGillespieNew.vtEpsilonW(sVoxel,2) = skJump_2;

end

% [x] OK: MK_correction: partially rewritten
function [rGillespieNew] =...
    G_Jump_in(sVoxel,rGillespie,sParticleOffset)
%
%   /// Description of input & output of the function /// 
%
%   MK_correction: revised the function.
%

%% Add the particle that just entered the voxel.
rGillespieNew = rGillespie;

if (rGillespieNew.sOffset_k > 0)

%(3.c.i) Update the cells Offset. 
%   Current time = rGillespie.st
%   time of entry = rGillespieNew.vOffset(sVoxel,2)
%   Offset decay rate constant = rGillespieNew.sOffset_k
rGillespieNew.vOffset(sVoxel,1) = rGillespieNew.vOffset(sVoxel,1) ...
    *exp( -rGillespieNew.sOffset_k * ...
    (rGillespieNew.st-rGillespieNew.vOffset(sVoxel,2) ));

else
    rGillespieNew.vOffset(sVoxel,1) = 0.0;
    sParticleOffset = 0.0;
end

%(3.c.ii) Add the incoming offset
rGillespieNew.vOffset(sVoxel,1) =  ...
    (rGillespieNew.mN(sVoxel,1)*1.0)/(rGillespieNew.mN(sVoxel,1)+1.0)*...
    rGillespieNew.vOffset(sVoxel,1)...
   +(1.0)/(rGillespieNew.mN(sVoxel,1)+1.0)*sParticleOffset;

%(3.c.iii) Reset the offset time to current time
rGillespieNew.vOffset(sVoxel,2) = rGillespie.st;

%(3.c.iv) Update Numbers
rGillespieNew.mN(sVoxel,1) = ...
    rGillespieNew.mN(sVoxel,1) + 1;

%% Recomputing The Jump Times:
% Calculating the effective sh_.
%    Positive Offset means smaller sh_, means jump to the right is more
%    likely.
rGillespieNew.sh_(sVoxel) = rGillespieNew.sh - ...
    2*rGillespieNew.vOffset(sVoxel,1);

% Cut-off such that h_>h_min>0.
if(rGillespieNew.sh_(sVoxel)<rGillespieNew.shMin)
    rGillespieNew.sh_(sVoxel) = rGillespieNew.shMin;
else
    if (rGillespieNew.sh_(sVoxel)>2-rGillespieNew.shMin)
        rGillespieNew.sh_(sVoxel) =2- rGillespieNew.shMin;
    end
end

% Local jump rate constant (if connected:
% to the left: 
skJump_1 = rGillespieNew.sD*rGillespieNew.mIcConnects(sVoxel,1)/...
    (2-rGillespieNew.sh_(sVoxel))^2;
% to the right:
skJump_2 = rGillespieNew.sD*rGillespieNew.mIcConnects(sVoxel,2)/...
    (rGillespieNew.sh_(sVoxel))^2;
% As noted above, due to division by sh_, for sh_ < 0.5 jump to the right
% is more likely.

%% Redraw the putative (jump) times.
rGillespieNew.vtEpsilon(sVoxel,1) = rGillespieNew.st+1.0/...
    (rGillespieNew.mN(sVoxel,1)*(skJump_1+skJump_2)).*exprnd(1,1);

% And store the jump weights:
rGillespieNew.vtEpsilonW(sVoxel,1) = skJump_1;
rGillespieNew.vtEpsilonW(sVoxel,2) = skJump_2;
end

% [x] OK: MK_correction: some corrections
function [rDataNew] = P_Create_particle(sVoxel,sDir,sParticleOffset,rData)
%
%   /// Description of input & output of the function /// 
%
%   MK_correction: revised the function.
%

%% Position of the transfered particle in the OmegaM domain and set active.
rDataNew = rData;

% Compute the resulting h_:
% (note that h=1.0 here)
sh_ = 1.0 - 2.0*sParticleOffset;

[~,sIx] = min(rDataNew.vActiveP);
rDataNew.vActiveP(sIx) = boolean(1);

if(sDir==2)
    rDataNew.vX(sIx) = sVoxel + 1 + ...
        min(sh_,1)*rand(1,1);
elseif(sDir==1)
    rDataNew.vX(sIx) = sVoxel - min(sh_,1)*rand(1,1);
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