%% Hybrid model (with offset)
%
%   Test a non-uniform particle distribution
%
%   This test uses a starting delta distribution for the particles, 
%   simulates the imroved method several times, and plots the normalized
%   distribution compared to the theoretical (expected) one.
%

%% Definitions for the experiment.
% Number of particles.
sNParticles = 10000;
% Starting position (delta distribution).
x0 = 21.0;
% Initialize particles position with a delta distribution and the species.
vX = ones(sNParticles,1).*x0;
vM = ones(sNParticles,1);
sNM = 1;

% vActiveP denotes with a boolean if the particle is in OmegaM.
vActiveP = boolean(ones(sNParticles,1));
% ["forbidden region" and OmegaC,OmegaC and OmegaM,OmegaM and "forbidden
% region"]
vBoundaries = [1.0,18.0,41.0];
% Duration of the simulation
st = 0; stEnd = 30; st0 = 0.0; sDeltaT = 0.01;
vT = st:sDeltaT:(stEnd-st); vTimes = [st,stEnd,st0,sDeltaT,vT];
% Number of trials the experiment is run.
sNTrials = 10;
% Bin width.
sbinWidth = 1.0;
% Number of voxels
sNK = vBoundaries(3)-vBoundaries(1);

% Offset decay in each voxel
OffsetDecayTunings = [0.0,8.5];

%% Simulations
for jj=1:length(OffsetDecayTunings)
mVoxelN = zeros(sNK,sNTrials);
mMeanPosition = zeros((1/sDeltaT)*stEnd,sNTrials);

tic;
for i = 1:sNTrials
    
    [rMytest2,vcenters,vVoxelN,vMeanPosition] = ...
        myTRM_1D_withOffset...
        (vX,vActiveP,vBoundaries,vTimes,OffsetDecayTunings(jj));
    
    if(sNTrials>1)
        mVoxelN(:,i) = vVoxelN;
        mMeanPosition(:,i) = vMeanPosition(:,1);
    end
end

sGillespieRuntime = toc;
fprintf('\nHybrid algorithm runtime %f s \n',sGillespieRuntime);
fprintf('\n')

if(sNTrials>1)
    vMeanVoxelN = mean(mVoxelN')';
    vMeanMeanPosition = mean(mMeanPosition')';
end

if(jj==1)
    rA = rMytest2;
else
    rB = rMytest2;
end
end

%% Compare to the normal distribution and expected distribution
if(1)    
    %vXTemp = rMytest2.vXNew(find(rMytest2.vActivePNew(:)==boolean(1)));
    %mNTemp = rMytest2.mmN(find(rMytest2.vvIcisGillespie==boolean(1)),end);
    
        
    vXTempA = rA.vXNew(find(rA.vActivePNew(:)==boolean(1)));
    mNTempA = rA.mmN(find(rA.vvIcisGillespie==boolean(1)),end);
        
    vXTempB = rB.vXNew(find(rB.vActivePNew(:)==boolean(1)));
    mNTempB = rB.mmN(find(rB.vvIcisGillespie==boolean(1)),end);
    
    % Generate the positions of the Gillespie particles.
    vPosGillA = [];
    for jj = 1:length(mNTempA)
       vPosGillA = [vPosGillA;ones(mNTempA(jj),1).*vcenters(jj) + ...
           rA.vOffset(jj)]; 
    end
    
    vPosGillB = [];
    for jj = 1:length(mNTempB)
       vPosGillB = [vPosGillB;ones(mNTempB(jj),1).*vcenters(jj) + ...
           rB.vOffset(jj)]; 
    end 
    
    [px_particleA,XA] = hist([vPosGillA;vXTempA],vcenters);
    [px_particleB,XB] = hist([vPosGillB;vXTempB],vcenters);
    
    
    % Normalize the curve such that the (numerical) integral over the whole
    % domain yields 1.0.
    px_particleA(:) = px_particleA(:)./sNParticles/sbinWidth;
    px_particleB(:) = px_particleB(:)./sNParticles/sbinWidth;
    
    sD = 1;
    px_theory = zeros(length(XA),1);
    for i = 1:length(XA)
        px_theory(i)=1/sqrt(4*pi*sD*stEnd)*...
            exp(-0.5*((XA(i)-x0)/sqrt(2*sD*stEnd))^2);
    end
    
    % Plot the particle - and the theoretical distributions.
    figure(2); clf; 
    plot(XA,px_particleA,'b','LineWidth',2); hold on; 
    plot(XB,px_particleB,'k','LineWidth',2); hold on; 
    plot(XA,px_theory,'r','LineWidth',2); grid on; 
    xlim([XA(1),XA(end)]);
    
    hXlabel = xlabel('x'); hYlabel = ylabel('Probability (normalized)');
    hLegend = legend('Obtained distribution',...
        'Obtained distribution (no offset)',...
        'Theoretical (expected) distribution'); 
    hTitle = title('Comparison distributions');
    
    set([hXlabel, hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
    
    hold off;
end
