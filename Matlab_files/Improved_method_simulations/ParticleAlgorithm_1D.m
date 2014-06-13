function [rDataN,rGillespieN] = ParticleAlgorithm_1D(rData,rGillespie)
%
%   Particle Algorithm (Brownian motion) (1D)
%
%   This function recalculates the positions of the active particles in the
%   Particle-domain. They move according to the Brownian motion.
%   Important to note are the conditions for particles abandoning OmegaM
%   and entering OmegaC, the Gillespie domain.
%

%% Constants and definitions

% Fix sDeltaX
sDeltaX = sqrt(2*rGillespie.sD*rGillespie.sDeltaT); 

%% Calculate new positions according to the Brownian equation.
vXOld = rData.vX;
vXNew = zeros(rData.sNParticles,1);
for j = 1:rData.sNParticles
    
    % If particle in the OmegaM domain.
    if(rData.vActiveP(j))
        
        vXNew(j) = vXOld(j) + sDeltaX*randn(1,1);
        
        vP = [vXNew(j),floor(vXNew(j))];
        
        % Check boundary condition between OmegaC and forbidden region 1
        % and OmegaM and forbidden region 2.
        % If out of bounds, particle jumps back to previous position.
        if (vP(1)>rGillespie.sxBoundaryFR2||vP(1)<rGillespie.sxBoundaryFR1)
            vP(1) = vXOld(j);
            vP(2) = floor(vXOld(j));
            
            % Store the obtained value.
            rData.vX(j) = vP(1);
  
        % Particle jumps from OmegaM to OmegaC
        elseif(rGillespie.vIcisGillespie(vP(2)))
            % Destroy particle in OmegaM domain.
            [rDataNew,sParticleOffset] =...
                P_Destroy_particle(j,rData,vP);
            rData = rDataNew;
                        
            % Create particle in OmegaC domain.
            [rGillespieNew] =...
                G_Jump_in(vP(2),rGillespie,sParticleOffset);
            rGillespie = rGillespieNew;
        
        % Store the obtained value because particle remains in OmegaM.
        else
            rData.vX(j) = vP(1);
        end
        
    % If particle in the OmegaC domain -> Do not move.
    else
        rData.vX(j) = vXOld(j); 
    end
    
end

rDataN = rData;
rGillespieN = rGillespie;

end

% [x] OK
function [rDataNew,sParticleOffset] = ...
    P_Destroy_particle(j,rData,vP)
%
%   /// Description of input & output of the function /// 
%

%% Delete the particle from the OmegaM domain.
rDataNew = rData;
rDataNew.vActiveP(j) = boolean(0);

%% Particle has a perfectly known position at the jump time, thus offset
%  can directly be calculated.

% Also have to allow negative offsets.
sParticleOffset = vP(1) - (vP(2)+0.5); 
% ==> Offset is positive if particle is on the right of the center of the 
%     voxel.

end

% [x] OK
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
    *exp( -rGillespieNew.sOffset_k * (rGillespieNew.st-rGillespieNew.vOffset(sVoxel,2) ));

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
fprintf('cut-off rGillespieNew.sh_(sVoxel) = %f \n',rGillespieNew.sh_(sVoxel));

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