%% Particle model
function [vvX] = offsetTest_1D(sPlots,sSimulation);
%
%   sPlots       : boolean, 1 outputs the plots
%   sSimulation  : boolean, 1 pauses some plots for simulation
%
%   vvX          : double matrix, contains x-values for the different   
%                  starting x0's (columnwise).

%% Constants and definitions
% Number of particles.
sNParticles = 10000;
% Starting positions (delta-distribution).
x0s = 0.5:0.1:0.9; 
% Here vActiveP irrelevant.
vActiveP = boolean(ones(sNParticles,1));
% Define mean and variance vectors.
sMean = zeros(2,1); sVar = zeros(2,1);
% Boundaries to the forbidden regions.
sxBoundaryFR1 = 0.0; sxBoundaryFR2 = 1.0;
% Center of the voxel.
sCenter = (sxBoundaryFR2 - sxBoundaryFR1)/2;
% Diffusion coefficient D = 1 [length_unit^2/time_unit].
sD = 1.0;
% Duration of the simulation / Time step / Time vector.
st = 0; stEnd = 1.0; sDeltaT = 0.01; vT = st:sDeltaT:(stEnd-st);
% Decaying constant for the exponentials.
c_2 = 8.5;

%% Simulations
bins = 0:0.01:1;
vvX = zeros(sNParticles,length(x0s));

for j = 1:length(x0s)
    
    if(sPlots==1 && sSimulation==1)
        figure((j-1)*2+1); clf;
    end
    
    x0 = x0s(j);
    vX = ones(sNParticles,1).*x0;
    sMean(j,1) = mean(vX); sVar(j,1) = var(vX);
    
    rData = struct('vX',vX,'vActiveP',vActiveP,'sNParticles',sNParticles);
    rGillespie = struct('sxBoundaryFR1',sxBoundaryFR1,'sxBoundaryFR2',...
        sxBoundaryFR2,'sDeltaT',sDeltaT,'sD',sD,'st',st,'stEnd',stEnd);
    
    %% Simulation
    ix = 1;
    vMeanOffset = zeros(length(vT),1);
    
    vOffset = rData.vX - sCenter.*ones(length(rData.vX),1);
    vMeanOffset(ix) = sum(vOffset)/length(vOffset);
    
    while(rGillespie.st<rGillespie.stEnd)
        
        % Plot the particle concentrations as histogram plot
        if(sPlots==1 && sSimulation==1)
            hist(rData.vX,bins); grid on; xlim([0,1]);
            xlabel('x'); ylabel('Number of particles');
            title('Domain is [0,1]');
            
            pause(0.3);
        end
        
        [vXNew,rGillespieNew] = particle1D(rData,rGillespie);
        
        rData.vX = vXNew;
        rGillespie = rGillespieNew;
        
        ix = ix + 1;
        
        vOffset = rData.vX - sCenter.*ones(length(rData.vX),1);
        vMeanOffset(ix) = sum(vOffset)/length(vOffset);
    end
   
    vvX(:,j) = rData.vX;
    
    sMean(j) = mean(vvX(:,j));
    sVar(j) = var(vvX(:,j));
    
    % Plot the exponential fitting functions
    if(sPlots==1)
        
        expFit(j,vT,vMeanOffset,x0,sCenter,c_2);
        
    end
    
end

fprintf('\n\n');
disp('The fitting function has the form:');
disp('f(t) = c_1*exp(-c_2*t)');
fprintf('\n');
disp('Note that c_1 = x0 - voxel_center.');
disp('The constant for the exponentials that seems to work for all the');
disp(['cases (except the first one, x_0 = 0.5) is c_2 = ',num2str(c_2)]);
fprintf('\n\n');

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

function expFit(j,vT,vMeanOffset,x0,sCenter,c_2)
figure(j*2); clf;

plot(vT,vMeanOffset,'ro','LineWidth',3); grid on; hold on;

if(j*2==2)
    figure(j*2); hold on;
    plot(vT,(x0-sCenter)*exp(-c_2*vT),'k','LineWidth',3); hold off;
    hXlabel = xlabel('time t'); 
    hYlabel = ylabel('offset with respect to x_0');
    hTitle = title(['Evolution of the offset for x_0 = ',num2str(x0)]);
    hLegend = legend('Data measurements',...
        'Fit function {0.0{\ite}^{-8.5*t}}');
    
    set([hXlabel,hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
end

if(j*2==4)
    figure(j*2); hold on;
    plot(vT,(x0-sCenter)*exp(-c_2*vT),'k','LineWidth',3); hold off;
    hXlabel = xlabel('time t');
    hYlabel = ylabel('offset with respect to x_0');
    hTitle = title(['Evolution of the offset for x_0 = ',num2str(x0)]);
    hLegend = legend('Data measurements',...
        'Fit function {0.1{\ite}^{-8.5*t}}');
    
    set([hXlabel,hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
end

if(j*2==6)
    figure(j*2); hold on;
    plot(vT,(x0-sCenter)*exp(-c_2*vT),'k','LineWidth',3); hold off;
    hXlabel = xlabel('time t');
    hYlabel = ylabel('offset with respect to x_0');
    hTitle = title(['Evolution of the offset for x_0 = ',num2str(x0)]);
    hLegend = legend('Data measurements',...
        'Fit function {0.2{\ite}^{-8.5*t}}');
    
    set([hXlabel,hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
end

if(j*2==8)
    figure(j*2); hold on;
    plot(vT,(x0-sCenter)*exp(-c_2*vT),'k','LineWidth',3); hold off;
    hXlabel = xlabel('time t');
    hYlabel = ylabel('offset with respect to x_0');
    hTitle = title(['Evolution of the offset for x_0 = ',num2str(x0)]);
    hLegend = legend('Data measurements',...
        'Fit function {0.3{\ite}^{-8.5*t}}');
    
    set([hXlabel,hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
end

if(j*2==10)
    figure(j*2); hold on;
    plot(vT,(x0-sCenter)*exp(-c_2*vT),'k','LineWidth',3); hold off;
    hXlabel = xlabel('time t');
    hYlabel = ylabel('offset with respect to x_0');
    hTitle = title(['Evolution of the offset for x_0 = ',num2str(x0)]);
    hLegend = legend('Data measurements',...
        'Fit function {0.4{\ite}^{-8.5*t}}');
    
    set([hXlabel,hYlabel,hLegend],'FontSize',16);
    set(hTitle,'FontSize',18,'FontWeight','bold');
end

end