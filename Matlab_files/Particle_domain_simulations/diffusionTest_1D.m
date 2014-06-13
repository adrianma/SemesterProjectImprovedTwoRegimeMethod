%% Particle model
function [vvX] = diffusionTest_1D(sPlots,sSimulation);
%
%   sPlots       : boolean, 1 outputs the plots.
%   sSimulation  : boolean, 1 pauses some plots for simulation.
%
%   vvX          : double vector, contains x-values for the different
%                  final positions.
%

%% Constants and definitions
% Number of particles.
sNParticles = 10000;
% Starting position (delta distribution).
x0 = 0.0;
% Boundaries to the forbidden regions.
sxBoundaryFR1 = -1000.0; sxBoundaryFR2 = 1000.0;
% Diffusion coefficient sD = 1 [length_unit^2/time_unit]
sD = 1.0;
% Duration of the simulation / Time step / Time vector.
st = 0; stEnd = 10.0; sDeltaT = 0.01; vT = st:sDeltaT:(stEnd-st);
% Steps
sN = (stEnd-st)/sDeltaT;
% Number of trials
sNTrials = 10;

%% Simulations
sMu = zeros(length(vT),1); sVarObserved = zeros(length(vT),1);
sbinWidth = 0.25;
bins = -50:sbinWidth:50;

vX = ones(sNParticles,1).*x0;
tic;
rData = struct('vX',vX,'sNParticles',sNParticles);
rGillespie = struct('sxBoundaryFR1',sxBoundaryFR1,'sxBoundaryFR2',...
    sxBoundaryFR2,'sDeltaT',sDeltaT,'sD',sD,'st',st,'stEnd',stEnd);

%% Simulation
ssMu = zeros(length(sNTrials),1);
ssVarObserved = zeros(length(sNTrials),1);

tic;
for jj = 1:sNTrials
    while(rGillespie.st<=rGillespie.stEnd)
        
        % Plot the particle concentrations as histogram plot
        if(sPlots==1 && sSimulation==1)
            figure(1);
            hist(rData.vX,bins); grid on;
            
            pause(0.3);
        end
        
        % Particle algorithm.
        [vXNew,rGillespieNew] = particle1D(rData,rGillespie);
        
        rData.vX = vXNew;
        rGillespie = rGillespieNew;
        
    end
    
    [sMu,sVarObserved] = normfit(rData.vX);
    sVarObserved = sVarObserved^2;
    sGillespieRuntime = toc;
    fprintf('\nHybrid algorithm runtime %f s \n',sGillespieRuntime);
    fprintf('\n');
    
    ssMu(jj) = sMu;
    ssVarObserved(jj) = sVarObserved;
    
end

sParticleRuntime = toc;
if(sSimulation==0)
    fprintf('\n \n');
    fprintf('Particle runtime %f s \n\n',sum(sParticleRuntime));
end

% Compute the mean for al sNTrials.
ssMu_mean = mean(ssMu);
ssVarObserved_mean = mean(ssVarObserved);

disp('HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH');
disp(['Mean of the mean for all simulations: ',num2str(ssMu_mean)]);
disp(['Mean of the observed variance for all simulations :',...
    num2str(ssVarObserved_mean)]);
disp(['Mean of the observed std for all simulations :',...
    num2str(sqrt(ssVarObserved_mean))]);
disp('HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH');
fprintf('\n');

% Plot the particle concentration as histogram plot (last iteration).
if(sPlots==1 && sSimulation==1)
    hist(rData.vX,bins); grid on;
    
    hXlabel = xlabel('x'); hYlabel = ylabel('Number of particles');
    hLegend = legend('Obtained distribution',...
        'Theoretical (expected) distribution'); 
    hTitle = title('Experiment for 10000 particles');
    
    set([hXlabel, hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
    
    pause(0.3)
end

% Store output
vvX = rData.vX;

%% Compute the (effective) diffusion coefficiont sDEff
sDEff = sVarObserved(end)/(2*rGillespie.sDeltaT*sN);
disp(['D_eff = ',num2str(sDEff),' compared to D = ',...
    num2str(rGillespie.sD)]);
fprintf('\n \n');

%% Compare to the normal distribution and expected distribution
if(1)
    % Declare vector with mean sMu and variance vvNormal.
    if(1)
        vvNormal = sMu(end) + sVarObserved(end).*randn(sNParticles,1);
        normalPlot(sPlots,rData.vX,vvNormal,bins);
    end
    
    [px_particle,X] = hist(rData.vX,bins);
    
    % Normalize the curve such that the (numerical) integral over the whole
    % domain yields 1.0.
    px_particle(:) = px_particle(:)./sNParticles/sbinWidth;
    
    px_theory = zeros(length(X),1);
    for i = 1:length(X)
        px_theory(i)=1/sqrt(4*pi*rGillespie.sD*rGillespie.stEnd)*...
            exp(-0.5*((X(i)-x0)/sqrt(2*rGillespie.sD*rGillespie.stEnd))^2);
    end
    
    % Plot the particle - and the theoretical distributions.
    figure(2); clf; 
    plot(X,px_particle,'LineWidth',2); hold on; 
    plot(X,px_theory,'r','LineWidth',2); grid on; 
    
    hXlabel = xlabel('x'); hYlabel = ylabel('Probability (normalized)');
    hLegend = legend('Obtained distribution',...
        'Theoretical (expected) distribution'); 
    hTitle = title('Comparison distributions');
    
    set([hXlabel, hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
        
    I = sum(sbinWidth.*px_particle);
    disp(['The value of the integral is I = ',num2str(I),',as expected.']);
    
end
end

function [vXNew,rGillespieNew] = particle1D(rData,rGillespie)

% Fix sDeltaX
sDeltaX = sqrt(2*rGillespie.sD*rGillespie.sDeltaT);

%% Calculate new positions according to the Brownian equation.
vXOld = rData.vX;
vXNew = zeros(rData.sNParticles,1);
for j = 1:rData.sNParticles
    
    vXNew(j) = vXOld(j) + sDeltaX*randn(1,1);
    
    %   vP = [x_position,chemical species,voxel_position]
    vP = [vXNew(j),floor(vXNew(j))];
    
    % Check boundary condition between OmegaM and forbidden region 1
    % and OmegaM and forbidden region 2.
    % If out of bounds, particle jumps back to previous position.
    if (vP(1)>rGillespie.sxBoundaryFR2||vP(1)<rGillespie.sxBoundaryFR1)
        vP(1) = vXOld(j);
        vP(2) = floor(vXOld(j));
        
        % Store the obtained value.
        vXNew(j) = vP(1);
        
    else
        % Store the obtained value.
        vXNew(j) = vP(1);
    end
    
end

rGillespie.st = rGillespie.st + rGillespie.sDeltaT;
rGillespieNew = rGillespie;

end

function normalPlot(sPlots,vvX,vvNormal,bins)

% Plot the particles positions.
if(sPlots==1)
    figure(3); clf;
    %subplot(2,1,1);
    hist(vvX,bins); grid on; xlim([bins(1),bins(end)]);
    
    hXlabel = xlabel('x'); hYlabel = ylabel('Number of particles');
    hTitle = title('Experiment for 10000 particles');
    
    set([hXlabel, hYlabel],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
    
    
    % Plot the normal distribution.
%     figure(3);
%     subplot(2,1,2);
%     hist(vvNormal,bins); grid on; xlim([bins(1),bins(end)]);
%     xlabel('x'); ylabel('Number of particles');
%     title('Normal distribution (for comparison)');
end

end
