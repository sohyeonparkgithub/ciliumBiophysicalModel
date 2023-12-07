%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modelSimulation.m
% Last updated on 2023-09-22 by Sohyeon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;


dataName = [datestr(now)];
saveName = dataName;


modelParameter;

%% Parameters to be learned

% SMO-IFT complex on and off % value obtained from Milenkovic, PNAS, 2015
kOn = 4; %param1
kOff = 2; %param2
%sterol oxisterol binding
kBind = 3; %param3 % within 8 seconds period, only one binding/unbinding event happened.
kUnbind = 1; %param4

v0 = 0.4; %0.4-0.5 µm/sec %Milenkovic, PNAS, 2015 %param5

xiMemDefault = 0.0158; %param6
xiMemTZ = 10*xiMemDefault; %known value not available %param7
uMemTZ = 100*kBT; %known value not available %param8

loc_error = 10^-5; %param9

%%
D0 = kBT/xiMemDefault; %0.26µm^2/s %Diffusion coefficient
threshold = 2*v0*dt+sqrt(D0*dt); %Doi model

%% Initialization
xState = NaN(maxSmo,nTime);

nowIFT = 0; % now-- tracks how many players are around now
nowSmo = 0;

numIFT = 0; %num-- how many players are created so far  (total counts)
numSmo = 0;

IFT = NaN(maxIFT,nTime);
smo = NaN(maxSmo,nTime);

%% Compute transition probability

pTrans = zeros(nState,nState); %transition probability

%states: 1.adv 2.diff 3. confinement(bind)

%transition probability
pTrans(1,2) = 1 - exp(-kOff*dt);
pTrans(1,3) = 0;
pTrans(1,1) = 1 - (pTrans(1,2)+pTrans(1,3));

pTrans(2,1) = (kOn/(kOn+kBind))*(1 - exp(-(kOn+kBind)*dt));
pTrans(2,3) = (kBind/(kOn+kBind))*(1 - exp(-(kOn+kBind)*dt));
pTrans(2,2) = 1 - (pTrans(2,1)+pTrans(2,3));

pTrans(3,1) = 0;
pTrans(3,2) = 1 - exp(-kUnbind*dt);
pTrans(3,3) = 1 - (pTrans(3,1)+pTrans(3,2));

%% energy barrier, viscosity, diffusion
uMem = @(location) uMemTZ.*exp(-(location-tzMiddle).^2./(2*(tzHalfWidth).^2));
uMemDeriv = @(location) -(exp(-(location - tzMiddle).^2/(2.*tzHalfWidth.^2)).*(uMemTZ).*(2.*location - 2.*tzMiddle))./(2.*tzHalfWidth.^2);

xiMem = @(location) xiMemDefault-(xiMemDefault-xiMemTZ).*exp(-(location-tzMiddle).^2./(2*(tzHalfWidth).^2));
xiCyto = @(location) xiCytoDefault-(xiCytoDefault-xiCytoTZ).*exp(-(location-tzMiddle).^2./(2*(tzHalfWidth).^2));
uCyto = @(location) uCytoDefault-(uCytoDefault-uCytoTZ).*exp(-(location-tzMiddle).^2./(2*(tzHalfWidth).^2));
uCytoDeriv = @(location) -(exp(-(location - tzMiddle).^2./(2.*tzHalfWidth.^2)).*(uCytoTZ - uCytoDefault).*(2.*location - 2.*tzMiddle))./(2.*tzHalfWidth.^2);
uCyto_uMem = @(location) (uCytoDefault) - (uCytoDefault -(uCytoTZ + uMemTZ)).*exp(-(location-tzMiddle).^2./(2*(tzHalfWidth).^2));
uDeriv = @(location) -(exp(-(location - tzMiddle).^2./(2.*tzHalfWidth.^2)).*(2.*location - 2.*tzMiddle).*(uCytoTZ - uCytoDefault + uMemTZ))./(2.*tzHalfWidth.^2);

denom = @(location) (xiMem(location) + xiCyto(location) + numMotors .* fStall./v0);
DVal = @(location) kBT./xiMem(location);
vIFT = @(location) (numMotors .* fStall - uCytoDeriv(location))./(xiCyto(location) + numMotors .* fStall./v0);
vSMO = @(location) - uMemDeriv(location)./xiMem(location);
vComplex = @(location) (numMotors .* fStall - uDeriv(location)) ./ denom(location);

%%
for iTime = 1:nTime-1

    if mod(iTime,1000)==0
        disp(iTime)
        datetime
    end

    %% Creation of IFT and smo
    t = iTime * dt;
    probIFT = lambdaIFT * dt; %IFT creation probability
    probSmo = lambdaSmo * dt;

    if (nowIFT < maxIFT)
        if (rand < probIFT)
            if (numIFT < maxIFT)
                numIFT = numIFT + 1;
                nowIFT = nowIFT + 1;
                IFT(numIFT,iTime) = 0;
            end %end of num loop
        end %end of prob loop
    end %end of IFT creation if loop

    if (nowSmo < maxSmo)
        if (rand < probSmo)
            if (numSmo < maxSmo)
                numSmo = numSmo + 1;
                nowSmo = nowSmo + 1;
                smo(numSmo,iTime) = 0;
                xState(numSmo,iTime) = 2;
            end %end of num loop
        end %end of prob loop
    end %end of smo creation if loop

    %% Movement loop
    for iIFT = 1:numIFT
        if (~isnan(IFT(iIFT,iTime)))
            IFT(iIFT,iTime+1) = IFT(iIFT,iTime) + vIFT(IFT(iIFT,iTime)) * dt;
        end %end of IFT ~isnan
    end % end of IFT movement

    for iSmo = 1:numSmo
        if (~isnan(smo(iSmo,iTime)))
            % Dynamics
            if (xState(iSmo,iTime) == 2) %diffusion
                smo(iSmo,iTime+1) = smo(iSmo,iTime) + vSMO(smo(iSmo,iTime)) * dt + normrnd(0,sqrt(2*DVal(smo(iSmo,iTime))*dt));
            elseif (xState(iSmo,iTime) == 1) %advection
                smo(iSmo,iTime+1) = smo(iSmo,iTime) + vComplex(smo(iSmo,iTime)) * dt;
            elseif (xState(iSmo,iTime) == 3) %confinement
                smo(iSmo,iTime+1) = smo(iSmo,iTime);
            end % end of if state loop
            %% Poisson-driven switch
            draw = rand;
            for iIFT = 1:numIFT
                if (xState(iSmo,iTime) == 1) %advection to diffusion kOff
                    if draw < pTrans(1,2)
                        xState(iSmo,iTime+1) = 2;
                    else
                        xState(iSmo,iTime+1) = 1;
                    end
                end
            end

            if (xState(iSmo,iTime) == 2) %diffusion to confinement kBind
                if draw < pTrans(2,3)
                    xState(iSmo,iTime+1) = 3;
                else
                    xState(iSmo,iTime+1) = 2;
                end
            elseif (xState(iSmo,iTime) == 3) %confinement to diffusion kUnbind
                if draw < pTrans(3,2)
                    xState(iSmo,iTime+1) = 2;
                else
                    xState(iSmo,iTime+1) = 3;
                end
            end

            %% Collision based switch
            if (xState(iSmo,iTime) == 2)
                for iIFT = 1:numIFT
                    if (abs(smo(iSmo,iTime+1)-IFT(iIFT,iTime+1)) < threshold)
                        xState(iSmo,iTime+1) = 1;
                        smo(iSmo,iTime+1) = IFT(iIFT,iTime+1);
                        IFT(iIFT,iTime+1) = NaN;
                        break;
                    else
                        if draw < pTrans(2,3) %edited on 07-09-2023
                            xState(iSmo,iTime+1) = 3;
                        else
                            xState(iSmo,iTime+1) = 2;
                        end
                    end
                end % end of multiple ift for loop
            end

        end % end of smo movement if
    end %end of smo movement

    %% Removal of players
    for iIFT = 1:numIFT
        if (IFT(iIFT,iTime+1) < 0)
            IFT(iIFT,iTime+1) = NaN;
            nowIFT = nowIFT - 1;
        elseif (IFT(iIFT,iTime+1) > ciliaLength)
            IFT(iIFT,iTime+1) = NaN;
            nowIFT = nowIFT - 1;
        end %end of if
    end %end of IFT removal for loop
    %
    for iSmo = 1:numSmo
        if smo(iSmo,iTime+1) < 0
            smo(iSmo,iTime+1) = NaN;
            nowSmo = nowSmo - 1;
        elseif (smo(iSmo,iTime+1) > ciliaLength)
            smo(iSmo,iTime+1) = NaN;
            nowSmo = nowSmo - 1;
        end %end of if
    end %end of smo removal for loop
end %end of time loop

% Add measurement error
tArray = dt*(1:numel(smo(1,:)));

for iIFT = 1:numIFT
    for iLength = 1:length(tArray)
        if ~isnan(IFT(iIFT,iLength))
            measureError = normrnd(0,loc_error,1,1);
            IFT(iIFT,iLength) = IFT(iIFT,iLength) + measureError;
        end
    end

end

for iSmo = 1:numSmo
    for iLength = 1:length(tArray)
        if ~isnan(smo(iSmo,iLength))
            measureError = normrnd(0,loc_error,1,1);
            smo(iSmo,iLength) = smo(iSmo,iLength) + measureError;
        end
    end
end

%% ----- Plots -----------------
%Figure 1 shows kymograph for
%1. IFT (red), 2. smo(green), and IFT-smo complex(yellow)
figure(100), clf;
hold all
box on
% IFT
for iIFT = 1:numIFT
    plot(tArray,IFT(iIFT,:),'r-','LineWidth',2)
end

% smo
for iSmo = 1:numSmo
    plot(tArray,smo(iSmo,:),'g-','LineWidth',2)
end

ylim([0 ciliaLength])
xlim([0 tArray(end)])

rectangle('Position', [0, tzMiddle-tzHalfWidth, tArray(end), 2*tzHalfWidth], ...
    'FaceColor', [0.9290, 0.6940, 0.1250, 0.3],'EdgeColor', [0.9290, 0.6940, 0.1250, 0.3]);

%interval = 100;
%set(gca,'XTick',0:interval:nTime)
%set(gca,'XTickLabel',0:5:dt*nTime)
xlabel('Time [seconds]')
ylabel('Distance from the base to the tip [um]')
title('kymograph')
set(gca,'Color','w')
set(gcf,'position',[114 504 861 175])

%%
figure(101), clf;
hold all
box on

% IFT
for iIFT = 1:numIFT
    plot(tArray,IFT(iIFT,:),'r-','LineWidth',2)
end

% smo
for iSmo = 1:numSmo
    plot(tArray,smo(iSmo,:),'g-','LineWidth',2)
end

ylim([0 ciliaLength])
xlim([0 tArray(end)])

% rectangle('Position', [0, tzMiddle-tzHalfWidth, tArray(end), 2*tzHalfWidth], ...
%     'FaceColor', [0.9290, 0.6940, 0.1250, 0.3],'EdgeColor', [0.9290, 0.6940, 0.1250, 0.3]);

%interval = 100;
%set(gca,'XTick',0:interval:nTime)
%set(gca,'XTickLabel',0:5:dt*nTime)
xlabel('Time [seconds]')
ylabel('Distance from the base to the tip [um]')
title('kymograph')
set(gca,'Color','w')
set(gcf,'position',[114 504 861 175])

%% convert into a format to feed in mcmc code
[rows, ~] = find(~isnan(smo));
rowVec = unique(rows);
maxTrack = length(rowVec);

x = NaN(maxTrack,nTime);
obs_all = NaN(maxTrack,nTime);
xState_all = NaN(maxTrack,nTime);

for iTrack = 1:maxTrack
    obs_all(iTrack,:) = smo(rowVec(iTrack),:);
    xState_all(iTrack,:) = xState(rowVec(iTrack),:);
end
%% convert into a format to feed in mcmc code
obs_tmp2 = NaN(maxTrack,nTime);
xStateTmpVec = NaN(maxTrack,nTime);
for iTrack = 1:maxTrack
    obs_single = obs_all(iTrack,:);
    obs_tmp = obs_single(~isnan(obs_single));

    xState_single = xState_all(iTrack,:);
    xState_tmp = xState_single(~isnan(xState_single));

    len_obs_tmp = length(obs_tmp);

    if len_obs_tmp > 50
        obs_tmp2(iTrack,1:len_obs_tmp) = obs_tmp;
        xStateTmpVec(iTrack,1:len_obs_tmp) = xState_tmp(1:len_obs_tmp);
    end
end

[rows, columns] = find(~isnan(obs_tmp2));
rowVec = unique(rows);
maxTrack = length(rowVec);
dObs_tmp = zeros(maxTrack,nTime-1);

for iTrack = 1:maxTrack
    obs(iTrack,:) = obs_tmp2(rowVec(iTrack),:);
    xStateGT(iTrack,:) = xStateTmpVec(rowVec(iTrack),:);
    dObs_tmp(iTrack,1:nTime-1) = diff(obs(iTrack,1:nTime));
end


%% Parameter ground truth
GT = {'kOn','kOff','kBind','kUnbind','v0','xiMemDefault','xiMemTZ','uMemTZ','localization error'};
GTVal = [kOn; kOff; kBind; kUnbind; v0; xiMemDefault; xiMemTZ; uMemTZ; loc_error];

%% working directory (change here before running the code)
if ~exist('dataName','var')
    dataName = 'test';
end

cd('/Users/sohyeonpark/ciliumTransport/Results/biophysicalSimulationLearning')
% currDate = strrep(datestr(datetime), ':', '_');
% mkdir([currDate '_' dataName])

mkdir(dataName)
cd(['/Users/sohyeonpark/ciliumTransport/Results/biophysicalSimulationLearning/' dataName])

%% Write table of ground truth parameter values
% T = table(kOn, kOff, kBind, kUnbind, v0, xiMemDefault, xiMemTZ, uMemTZ, loc_error, 'VariableNames', {'kOn','kOff','kBind','kUnbind','v0','xiMemDefault','xiMemTZ','uMemTZ','localization error'})
T = table({'kOn';'kOff';'kBind';'kUnbind';'v0';'xiMemDefault';'xiMemTZ';'uMemTZ';'localizationError';},...
    [kOn;kOff;kBind;kUnbind;v0;xiMemDefault;xiMemTZ;uMemTZ;loc_error;])
writetable(T,'parameter.txt')


%% State sequence
figure(111);clf;
hold all
for iTrack = 1:maxTrack
    subplot(ceil(maxTrack/2),2,iTrack)
    plot(xStateGT(iTrack,:),'b-','LineWidth',2)
    xlabel('Time')
    ylabel('states')
    nonNanInd = find(~isnan(xStateGT(iTrack,:)));
    xlim([0 nonNanInd(end)])
    yticks([1 2 3])
    title(['track',num2str(iTrack),', Ground truth'])
end % multiple track for loop

%% Simulated kymograph, no rectangle on TZ
%1. IFT (red), 2. smo(green), and IFT-smo complex(yellow)
figure(101), clf;
hold all
box on
% IFT
% for iIFT = 1:numIFT
%     plot(tArray,IFT(iIFT,:),'r-','LineWidth',2)
% end

% smo
for iSmo = 1:numSmo
    plot(tArray,smo(iSmo,:),'g-','LineWidth',2)
end

ylim([0 ciliaLength])
xlim([0 tArray(end)])


% rectangle('Position', [0, tzMiddle-tzHalfWidth, tArray(end), 2*tzHalfWidth], ...
%     'FaceColor', [0.9290, 0.6940, 0.1250, 0.3],'EdgeColor', [0.9290, 0.6940, 0.1250, 0.3]);

%interval = 100; 
%set(gca,'XTick',0:interval:nTime)
%set(gca,'XTickLabel',0:5:dt*nTime)
xlabel('Time [seconds]')
ylabel('Distance from the base to the tip [um]')
title('kymograph')
% set(gca,'Color','k')
% set(gca,'Color','w')
set(gcf,'position',[114 504 861 175])
saveas(gcf,'BPsimulation','epsc')


%% Simulated kymograph
%1. IFT (red), 2. smo(green), and IFT-smo complex(yellow)
figure(100), clf;
hold all
box on
%IFT
for iIFT = 1:numIFT
    plot(tArray,IFT(iIFT,:),'r-','LineWidth',2)
end

% smo
for iSmo = 1:numSmo
    plot(tArray,smo(iSmo,:),'g-','LineWidth',2)
end

ylim([0 ciliaLength])
xlim([0 tArray(end)])


rectangle('Position', [0, tzMiddle-tzHalfWidth, tArray(end), 2*tzHalfWidth], ...
    'FaceColor', [0.9290, 0.6940, 0.1250, 0.3],'EdgeColor', [0.9290, 0.6940, 0.1250, 0.3]);

%interval = 100;
%set(gca,'XTick',0:interval:nTime)
%set(gca,'XTickLabel',0:5:dt*nTime)
xlabel('Time [seconds]')
ylabel('Distance from the base to the tip [um]')
title('kymograph')
% set(gca,'Color','k')
% set(gca,'Color','w')
set(gcf,'position',[114 504 861 175])
saveas(gcf,'BPsimulationTZ','epsc')



%% Clean up
clearvars -except obs dObs_tmp maxTrack GT GTVal xStateGT dataName saveName
save([saveName,'.mat'])