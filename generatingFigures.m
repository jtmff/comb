clear
%% Loading data
dataCa140 = load('bcl140_pacing_Ca');
dataVm140 = load('bcl140_pacing_V');
dataCa80alternans = load('bcl85_alternans_Ca');
dataCa80alternansSevere = load('bcl85_alternansSevere_Ca');

mkdir imgs
%% Normal Ca signal
minimaCa140 = 20:140:1000;
minimaCa140V2 = 65:140:1000;
minimaCa140V3 = 110:140:1000;

figure(1); clf;
h(1) = plot(dataCa140.avgTrace, 'LineWidth', 1.25);
% comb 1
hold on
for iMin = 1:length(minimaCa140)
    plot([minimaCa140(iMin), minimaCa140(iMin)], [1060, 1130], 'k--', 'LineWidth', 1.5);
end
h(2) = plot([minimaCa140(1), 1000], [1130, 1130], 'k--', 'LineWidth', 1.5);
hold off

% comb  2
hold on
for iMin = 1:length(minimaCa140V2)
    plot([minimaCa140V2(iMin), minimaCa140V2(iMin)], [1070, 1140], 'Color',[128 128 128]/255, 'LineWidth', 1.5);
end
h(3) = plot([minimaCa140V2(1), 1000], [1140, 1140], 'Color',[128 128 128]/255, 'LineWidth', 1.5);
hold off

% comb  3
hold on
for iMin = 1:length(minimaCa140V3)
    plot([minimaCa140V3(iMin), minimaCa140V3(iMin)], [1080, 1150], ':', 'Color',[64 64 64]/255, 'LineWidth', 1.5);
end
h(4) = plot([minimaCa140V3(1), 1000], [1150, 1150], ':', 'Color',[64 64 64]/255, 'LineWidth', 1.5);
hold off

mean1 = mean(dataCa140.avgTrace(minimaCa140));
mean2 = mean(dataCa140.avgTrace(minimaCa140V2));
mean3 = mean(dataCa140.avgTrace(minimaCa140V3));

ylim([1000, 1160]);
legend(h, {'Signal', 'Comb 1 (mean = 1051)', 'Comb 2 (mean = 1112)', 'Comb 3 (mean = 1075)'}, 'Location', 'southeast');
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 16 12]);
saveas(gcf, 'imgs/bcl140_pacing_Ca_comb.png'); 

%% Normal Vm signal
minimaV140 = combGetMinima(dataVm140.avgTrace, 140);
apBoundaries = minimaV140 - 20; % assuming upstroke is max 20 ms.

figure(2); clf;
plot(2.02e4 - dataVm140.avgTrace, 'LineWidth', 1.25);
hold on
for iBoundary = 1:length(minimaV140)
    plot([apBoundaries(iBoundary), apBoundaries(iBoundary)], [0 1200], 'k--');
end
hold off

legend('Signal', 'Boundaries');
xlabel('Time (ms)');
ylabel('Mem. pot. dye fluorescence');
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 16 12]);
saveas(gcf, 'imgs/bcl140_pacing_V.png'); 

%% Linear and nonlinear baseline drift, stepwise function added, and gradual amplitude reduction
% linear and nonlinear
linearDriftCa = dataCa140.avgTrace + (1:length(dataCa140.avgTrace))'/20;
nonlinearDriftCa = dataCa140.avgTrace + ((1:length(dataCa140.avgTrace))'/20).^1.5 -((1:length(dataCa140.avgTrace))'/150).^3 ;
minimaCa140linearDrift = combGetMinima(linearDriftCa, 140);
minimaCa140nonlinearDrift = combGetMinima(nonlinearDriftCa, 140);

% stepwise drift
stepSignal = dataCa140.avgTrace + [zeros(290+50,1); 500*ones(427-50,1); zeros(283,1)];
minimaCa140stepSignal = combGetMinima(stepSignal, 140);

% reducing amplitude

signalReducingAmplitude = dataCa140.avgTrace' .* (0.9998:-4e-4:0.6) + (1:length(dataCa140.avgTrace))/2.5;
minimaCa140amplReduction = combGetMinima(signalReducingAmplitude, 140);

figure(3); clf
subplot(2,2,1);
plot(linearDriftCa, 'LineWidth', 1.25)
hold on
for iMin = 1:length(minimaCa140linearDrift)
    plot([minimaCa140linearDrift(iMin), minimaCa140linearDrift(iMin)], [1000 1180], 'k--');
end
hold off
ylim([1000, 1180]);
legend({'Signal', 'Boundaries'}, 'Location', 'southeast');
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');
title('Linear drift');


subplot(2,2,2);
plot(nonlinearDriftCa, 'LineWidth', 1.25)
hold on
for iMin = 1:length(minimaCa140nonlinearDrift)
    plot([minimaCa140nonlinearDrift(iMin), minimaCa140nonlinearDrift(iMin)], [1025 1250], 'k--');
end
hold off
ylim([1025, 1250]);
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');
title('Nonlinear drift')

subplot(2,2,3)
plot(stepSignal, 'LineWidth', 1.25)
hold on
for iMin = 1:length(minimaCa140stepSignal)
    plot([minimaCa140stepSignal(iMin), minimaCa140stepSignal(iMin)], [1025 1750], 'k--');
end
hold off
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');
title('Step function added');
ylim([1025, 1750]);

subplot(2,2,4);
plot(signalReducingAmplitude, 'LineWidth', 1.25)
hold on
for iMin = 1:length(minimaCa140amplReduction)
    plot([minimaCa140amplReduction(iMin), minimaCa140amplReduction(iMin)], [800 1250], 'k--');
end
hold off
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');
title('Amplitude reduction');
ylim([800, 1250]);

set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 24 24]);
saveas(gcf, 'imgs/bcl140_pacing_drift.png'); 


%% with white noise
lowNoise = dataCa140.avgTrace + randn(size(dataCa140.avgTrace)) * 10;
midNoise = dataCa140.avgTrace + randn(size(dataCa140.avgTrace)) * 25;
highNoise = dataCa140.avgTrace + randn(size(dataCa140.avgTrace)) * 50;

minimaLowNoise = combGetMinima(lowNoise, 140);
minimaMidNoise = combGetMinima(midNoise, 140);
minimaHighNoise = combGetMinima(highNoise, 140);

figure(4);
subplot(3,1,1)
plot(lowNoise, 'LineWidth', 1.25);
hold on
for iMin = 1:length(minimaLowNoise)
    plot([minimaLowNoise(iMin), minimaLowNoise(iMin)], [900 1300], 'k--');
end
hold off
ylim([900 1300]);
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');
title('Noise sd = 10');
legend('Signal', 'Boundaries');

subplot(3,1,2);
plot(midNoise, 'LineWidth', 1.25);
hold on
for iMin = 1:length(minimaMidNoise)
    plot([minimaMidNoise(iMin), minimaMidNoise(iMin)], [900 1300], 'k--');
end
hold off
ylim([900 1300]);
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');
title('Noise sd = 25');

subplot(3,1,3);
plot(highNoise, 'LineWidth', 1.25);
hold on
for iMin = 1:length(minimaHighNoise)
    plot([minimaHighNoise(iMin), minimaHighNoise(iMin)], [900 1300], 'k--');
end
hold off
ylim([900 1300]);
title('Noise sd = 50');
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');

set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 24 24]);
saveas(gcf, 'imgs/bcl140_gaussianNoise.png');

%% Salt-and-pepper noise
% Case 1 - decently behaved salt and pepper noise
indicesMinimum = randi(1000,1,20);
indicesMaximum = randi(1000,1,20);
saltPepper = dataCa140.avgTrace;
saltPepper(indicesMinimum) = 950;
saltPepper(indicesMaximum) = 1150;
minimaSaltPepper = combGetMinima(saltPepper, 140);

% Case 2 - only one point of massive amplitude.
saltPepperEvil = dataCa140.avgTrace;
saltPepperEvil(500) = 0;
minimaSaltPepperEvil = combGetMinima(saltPepperEvil, 140);

% Plotting both
figure(5);
subplot(2,1,1);
plot(saltPepper, 'LineWidth', 1.25);
hold on
for iMin = 1:length(minimaSaltPepper)
    plot([minimaSaltPepper(iMin), minimaSaltPepper(iMin)], [850 1300], 'k--');
end
hold off
ylim([850 1250]);
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');

subplot(2,1,2);
plot(saltPepperEvil, 'LineWidth', 1.25);
hold on
for iMin = 1:length(minimaSaltPepperEvil)
    plot([minimaSaltPepperEvil(iMin), minimaSaltPepperEvil(iMin)], [0 1250], 'k--');
end
hold off
ylim([0 1250]);
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');
legend({'Signal', 'Boundaries'}, 'Location','southeast');

set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 24 20]);
saveas(gcf, 'imgs/bcl140_saltPepperNoise.png');

%% Alternans
minimaAlternans = combGetMinima(dataCa80alternans.avgTrace, 85);
minimaAlternansSevere = combGetMinima(dataCa80alternansSevere.avgTrace, 85);
figure(6);
subplot(2,1,1);
plot(dataCa80alternans.avgTrace, 'LineWidth', 1.25);
hold on
for iMin = 1:length(minimaAlternans)
    plot([minimaAlternans(iMin), minimaAlternans(iMin)], [1130 1200], 'k--');
end
hold off
ylim([1130 1200]);
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');
legend('Signal', 'Boundaries');

subplot(2,1,2);
plot(dataCa80alternansSevere.avgTrace, 'LineWidth', 1.25);
hold on
for iMin = 1:length(minimaAlternansSevere)
    plot([minimaAlternansSevere(iMin), minimaAlternansSevere(iMin)], [1050 1180], 'k--');
end
hold off
ylim([1050 1180]);
xlabel('Time (ms)');
ylabel('Ca^{2+} dye fluorescence');



set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 24 20]);
saveas(gcf, 'imgs/bcl85_alternans.png');
