function [] = PlotDesignResults()

% Load data
load('phasesep_T0.10.mat')

% Setup figures
close all
M = size(ext,1);
N = size(x,1);
h_ext = figure; hold on
h_error = figure; hold on

% Grab logarithmically-spaced frames
N_plot = 100;
plot_ind = unique(round(exp(linspace(0,log(M-1),N_plot)')));
plot_ind = [plot_ind;M];
N_plot = length(plot_ind);
colors = parula(N_plot);

% Plot target extinction
figure(h_ext)
plot(omega, ext_t, 'r', 'LineWidth', 3)

% Loop through iterations and plot extinction
for i = 1:N_plot
    %plot(omega, ext(plot_ind(i),:), 'Color', colors(i,:), 'LineWidth', 1)
end

% Format extinction
figure(h_ext)
xlabel('\omega')
ylabel('\sigma')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
xlim([0, 1])
set(gca, 'XDir', 'reverse')

% Plot error
figure(h_error)
plot(1:M, E, 'o')

% Format error
figure(h_error)
xlabel('iteration')
ylabel('error')
set(gca, 'XScale', 'log', 'YScale', 'log')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])

end