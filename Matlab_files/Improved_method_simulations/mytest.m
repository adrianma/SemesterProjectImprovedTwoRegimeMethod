%% Hybrid model (with offset)
%
%   This simulations correspond to the ones to be found on the written
%   report, sections 4.2 and 5.
%
%

%% Definitions for the experiment.
% Number of particles.
sNParticles = 2000;
% Initialize particles position and species.
% // Random uniform distribution //
vX = 1+rand(sNParticles,1).*40.0;
vM = ones(sNParticles,1);
sNM = 1;

% vActiveP denotes with a boolean if the particle is in OmegaM.
vActiveP = boolean(ones(sNParticles,1));
% ["forbidden region" and OmegaC,OmegaC and OmegaM,OmegaM and "forbidden
% region"]
vBoundaries = [1.0,31.0,41.0];
% Duration of the simulation
st = 0; stEnd = 30.0; st0 = 0.0; sDeltaT = 0.01;vT = st:sDeltaT:(stEnd-st);
vTimes = [st,stEnd,st0,sDeltaT,vT];
% Number of trials the experiment is run.
sNTrials = 2;
% Bin width.
sbinWidth = 1.0;
% Number of voxels
sNK = vBoundaries(3)-vBoundaries(1);

% Offset decay in each voxel
% if 0: then it is run without offset
OffsetDecayTunings = [2.5,0.0,8.5]; 

fprintf('\n');
disp(['Domain is bounded between ',num2str(vBoundaries(1)),...
    '.0 and ',num2str(vBoundaries(3)),'.0,']);
disp(...
    ['and the interface between Gillespie and Particle domain is at '...
    ,num2str(vBoundaries(2)),'.0.']);
fprintf('\n');
disp(['Run the algorithm for ',num2str(stEnd),...
    ' s. The time step is ',num2str(sDeltaT),' s.']);
fprintf('\n');
disp(['The whole experiment is run ',num2str(sNTrials),' times.']);
fprintf('\n');

%% Simulations
% Run the improved algorithm for the three different offset jumping rate
% constants.

sGillespieRuntime = zeros(length(OffsetDecayTunings),1);
vvv = zeros((1/sDeltaT)*stEnd + 1,length(OffsetDecayTunings));
vvv_t = zeros((1/sDeltaT)*stEnd + 1,length(OffsetDecayTunings));


tic;
for jj=1:length(OffsetDecayTunings)
    
    OffsetDecayTuning = OffsetDecayTunings(jj);
    
    mVoxelN = zeros(sNK,sNTrials);
    mMeanPosition = zeros(stEnd*100 + 1,sNTrials);
    
    for i = 1:sNTrials
        
        [~,vcenters,vVoxelN,vMeanPosition] = ...
            myTRM_1D_withOffset...
            (vX,vActiveP,vBoundaries,vTimes,OffsetDecayTuning);
        if(sNTrials>1)
            mVoxelN(:,i) = vVoxelN;
            mMeanPosition(:,i) = vMeanPosition(:,1);
        end
    end
    
    % Average the results of the 10 runs.
    if(sNTrials>1)
        vMeanVoxelN = mean(mVoxelN')';
        vMeanMeanPosition = mean(mMeanPosition')';
    end
    
    %% Plots for comparison for sNTrials realizations
    if(1)
        
        if(sNTrials>1)
            vvVoxelN = vMeanVoxelN;
            vvMeanPosition = vMeanMeanPosition;
        else
            vvVoxelN = vVoxelN;
            vvMeanPosition = vMeanPosition;
        end
        
        % Plot the particle distribution
        figure(jj); clf;
        stem(vcenters,vvVoxelN,'.','LineWidth',3); grid on;
        hold on;
        % Plot the interface line.
        line([vBoundaries(2),vBoundaries(2)],[0.0,70],...
            'Color','k','LineWidth',3);
        hXlabel = xlabel('Voxel');
        hYlabel = ylabel('Number of particles in given voxel');
        hTitle = title('Molecule distribution in space');
        hLegend = legend(['Offset decay k=',num2str(OffsetDecayTuning)]);
        set([hXlabel,hYlabel,hLegend],'FontSize',16);
        set(hTitle,'FontSize',18,'FontWeight','bold');
        
        vvv(:,jj) = vvMeanPosition(:,1);
        vvv_t(:,jj) = vMeanPosition(:,2);
        
        % Plot the mean position over time, for offset discussion.
        figure(jj+length(OffsetDecayTunings)); clf;
        plot(vMeanPosition(:,2),vvMeanPosition(:,1));
        xlim([st,stEnd]); grid on;
        hXlabel = xlabel('t');
        hYlabel = ylabel('Mean position offset');
        hTitle = title('Molecule mean position over time');
        hLegend = legend(['Offset decay k=',num2str(OffsetDecayTuning)]);
        set([hXlabel,hYlabel,hLegend],'FontSize',16);
        set(hTitle,'FontSize',18,'FontWeight','bold');
        
        
    end
end

% Print out the running time for the simulation.
sGillespieRuntime(jj) = toc;
fprintf('\nHybrid algorithm runtime %f s \n',sGillespieRuntime(jj));
fprintf('\n');

%% Plot all mean offsets in one figure for comparison.
if(1)
    figure;
    hold on; grid on;
    plot(vvv_t(:,1),vvv(:,1),'b','LineWidth',2);
    plot(vvv_t(:,2),vvv(:,2),'k','LineWidth',2);
    plot(vvv_t(:,3),vvv(:,3),'r','LineWidth',2);
    xlim([vvv_t(1,1),vvv_t(end,1)]);
    
    hXlabel = xlabel('t');
    hYlabel = ylabel('Mean position offset');
    hTitle = title('Comparison molecule mean position over time');
    hLegend = legend('k=2.5','k=0.0','k=8.5');
    set([hXlabel,hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
end
