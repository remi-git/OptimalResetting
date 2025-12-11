% This code plots the results of the associated python code [it should therefore be run after].
% It produces three subplots, (1) the measured mean first passage time as a
% function of the protocol duration (2) the measured work as a function of
% the protocol duration and (3) the time-energy bound, i.e. the measured
% MFPT as a function of the measured work.
% It does so for a series of potential stiffness, which includes the case
% studied in the main text of the paper, but also cases with low stiffness,
% where scales are not well separated anymore.

clear all
clc
dt = 1e-5; % time-steps
kB = 1.3806e-23; % Boltzmann constant
T = 300; % room temperature
gamma = 2.8e-8; % Stokes viscous drag
D = kB * T / gamma; % diffusion coefficient
x0 = 0.1e-6; % target
OptimalRate = 36; % resetting rate
Mfpt_inst =  1./OptimalRate .* (exp(sqrt(OptimalRate/D) * x0) - 1); % instantaneous MFPT
z = x0 * sqrt(OptimalRate / D);
CorrectionSecondMoment = 1 - 1 / (1 - exp(z)) + z^2 / (4 * (1 - cosh(z * sqrt(1 - exp(-z))))); % taking into account the absorbing boundary condition
AllKappaValue = flip([1, 5, 25, 50, 100]);
NCaseKappaPlot = numel(AllKappaValue);
Cmap = customcolormap([0 0.6 1],...
    [[0.1690, 0.4040, 0.6350]; [0.9690, 0.7000, 0.0860] ...
    ; [0.8200, 0.3900, 0.3920]], NCaseKappaPlot);
Symbols = ['d', '^', 's', 'o', 'v', '>'];
GreyBlue = [0.1490, 0.3180, 0.3570];
MinTimeArray = flip(linspace(0.005, 0.5, NCaseKappaPlot));
MaxTimeArray = flip(linspace(5, 35, NCaseKappaPlot));

figure()
subplot(2,2,1)
for k = 1:NCaseKappaPlot % plot the measured MFPT versus protocol duration
    KappaValue = AllKappaValue(k);
    Stiffness = KappaValue * 1e-6;
    AllMeanFptSim = load("AllMeanfpt"  + KappaValue +  ".mat").AllMeanfpt;
    AllErrFptSim = load("AllErrfpt"  + KappaValue +  ".mat").AllErrfpt;
    AllMeanFptInstSim = load("AllMeanfptInst"  + KappaValue +  ".mat").AllMeanfptInst;
    Alltf = gamma / Stiffness .* logspace(log10(MinTimeArray(k)), log10(MaxTimeArray(k)), numel(AllMeanFptSim));
    p1 = errorbar(1e3*Alltf, AllMeanFptSim./AllMeanFptInstSim, ...
        AllErrFptSim./AllMeanFptInstSim, 'vertical',...
        Symbols(k), 'linewidth', 2, 'Color', Cmap(k,:));
    xscale log
    yscale log
    p1.CapSize = 0;
    p1.MarkerSize = 8;
    hold on
end
tfArray = gamma / Stiffness * linspace(3e-3, 8, 10000);
ConstrainedProtocol_Mfpt = Mfpt_inst .* (1 + OptimalRate * tfArray); % constrained duration
loglog(1e3*tfArray, ConstrainedProtocol_Mfpt./Mfpt_inst, 'color', GreyBlue, 'linewidth', 2)
axis tight
xticks([1e-1, 1e0, 1e1, 1e2])
xticklabels({'10^{-1}', '10^0', '10^1', '10^2'})
yticks([1, 10])
yticklabels({'10^0', '10^1'})
xlabel('$ t_{\rm f}  ~\rm{[ms]}$', 'interpreter', 'latex')
ylabel('$ \langle \mathcal{T} \rangle / \langle \mathcal{T}_{\rm inst} \rangle$', 'interpreter', 'latex')
set(gca,'FontSize', 20, 'fontname', 'times');
box on
subplot(2,2,2)
for k = 1:NCaseKappaPlot % plot the measured work versus protocol duration
    KappaValue = AllKappaValue(k);
    Stiffness = KappaValue * 1e-6;
    AllMeanWorkSim = load("AllMeanWork" + KappaValue + ".mat").AllMeanWork;
    AllMeanFptSim = load("AllMeanfpt"  + KappaValue +  ".mat").AllMeanfpt;
    AllErrWorkSim = load("AllErrWork"  + KappaValue +  ".mat").AllErrWork;
    Alltf = gamma / Stiffness .* logspace(log10(MinTimeArray(k)), log10(MaxTimeArray(k)), numel(AllMeanFptSim));
    p1 = errorbar(1e3*Alltf, AllMeanWorkSim./(kB*T), ...
        AllErrWorkSim./(kB*T), 'vertical',...
        Symbols(k), 'linewidth', 2, 'Color', Cmap(k,:));
    xscale log
    yscale log
    p1.CapSize = 0;
    p1.MarkerSize = 8;
    hold on
end
ArrayWorkOptimal = CorrectionSecondMoment * 2 * kB * T ./ (OptimalRate .* tfArray); % cost optimal protocol
loglog(1e3*tfArray, ArrayWorkOptimal./(kB*T), 'color', GreyBlue, 'linewidth', 2)
axis tight
xticks([1e-1, 1e0, 1e1, 1e2])
xticklabels({'10^{-1}', '10^0', '10^1', '10^2'})
yticks([1, 100])
yticklabels({'10^0', '10^2'})
xlabel('$ t_{\rm f}  ~\rm{[ms]}$', 'interpreter', 'latex')
ylabel('$ \langle \Delta W \rangle  ~\rm{[k_B T]}$', 'interpreter', 'latex')
set(gca,'FontSize', 20, 'fontname', 'times');
box on
subplot(2,2,[3,4])
for k = 1:NCaseKappaPlot % plot MFPT verus work (time-energy bound)
    KappaValue = AllKappaValue(k);
    Stiffness = KappaValue * 1e-6;
    AllMeanWorkSim = load("AllMeanWork" + KappaValue + ".mat").AllMeanWork;
    AllMeanFptSim = load("AllMeanfpt"  + KappaValue +  ".mat").AllMeanfpt;
    AllErrWorkSim = load("AllErrWork"  + KappaValue +  ".mat").AllErrWork;
    AllErrFptSim = load("AllErrfpt"  + KappaValue +  ".mat").AllErrfpt;
    AllMeanFptInstSim = load("AllMeanfptInst"  + KappaValue +  ".mat").AllMeanfptInst;
    p1 = plot(AllMeanWorkSim/(kB*T), AllMeanFptSim./AllMeanFptInstSim, ...
        Symbols(k), 'linewidth', 3, 'Color', Cmap(k,:));
    xscale log
    yscale log
    p1.MarkerSize = 7;
    hold on
end
tfArray = gamma / min(AllKappaValue*1e-6) * linspace(0.01, 80, 10000);
ConstrainedProtocol_Mfpt = Mfpt_inst .* (1 + OptimalRate * tfArray); % constrained duration
ArrayWorkOptimal = CorrectionSecondMoment * 2 * kB * T ./ (OptimalRate .* tfArray); % cost optimal protocol
loglog(ArrayWorkOptimal/(kB*T), ConstrainedProtocol_Mfpt./Mfpt_inst, 'color', GreyBlue, 'linewidth', 2)
area(ArrayWorkOptimal/ (kB * T), ConstrainedProtocol_Mfpt./Mfpt_inst, 'FaceColor', GreyBlue, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
legend('\kappa = 100 (\sigma \approx \zeta / 10)', '\kappa = 50', ...
    '\kappa = 25', '\kappa = 5', '\kappa = 1 (\sigma \approx \zeta)', 'Analytical')
yticks([1, 10])
yticklabels({'10^0', '10^1'})
axis tight
xlim([0.2, 2e2])
ylim([0.9, 12])
xlabel('$ \langle \Delta W \rangle  ~\rm{[k_B T]}$', 'interpreter', 'latex')
ylabel('$ \langle \mathcal{T} \rangle / \langle \mathcal{T}_{\rm inst} \rangle$', 'interpreter', 'latex')
set(gca,'FontSize', 20, 'fontname', 'times');
set(gcf, 'Color', 'w');
set(gcf, 'units', 'centimeters', 'position', [1 1 25 20]);
set(gca, 'fontname', 'times'); 
exportgraphics(gcf, "ResettingVariousKappa.png",'Resolution', 200)
