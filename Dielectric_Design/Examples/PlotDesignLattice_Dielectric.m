function [] = PlotDesignLattice_Dielectric()

% Load data
load('test_2.mat')

% Setup figures
close all
M = size(ext,1);
m = size(params,1);
h_ext = figure; hold on
h_error = figure; hold on
colors = parula(M);

% Plot target extinction
figure(h_ext)
plot(omega, ext_t, 'r', 'LineWidth', 3)

% Loop through iterations and plot extinction
for i = 1:M
    plot(omega, ext(i,:), 'Color', colors(i,:), 'LineWidth', 1)
end

% Format extinction
figure(h_ext)
xlabel('\omega')
ylabel('\sigma')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
xlim([0.2, 1])
set(gca, 'XDir', 'reverse')

% Plot error
figure(h_error)
plot(1:M, E, 'o')

% Format error
figure(h_error)
xlabel('iteration')
ylabel('error')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1], 'YScale', 'log')

colors = parula(m);

% Plot parameters
figure
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'none')

nexttile(1); hold on
for i = 1:m
    plot(1:M, squeeze(params(i,1,:)), 'Color', colors(i,:), 'LineWidth', 1)
end
ylabel('a')
set(gca,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[3,1,1])

nexttile(2); hold on
for i = 1:m
    plot(1:M, squeeze(params(i,2,:)), 'Color', colors(i,:), 'LineWidth', 1)
end
ylabel('\omega_p')
set(gca,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[3,1,1])

nexttile(3); hold on
for i = 1:m
    plot(1:M, squeeze(params(i,3,:)), 'Color', colors(i,:), 'LineWidth', 1)
end
ylabel('\gamma')
set(gca,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[3,1,1])

nexttile(4); hold on
for i = 1:m
    plot(1:M, squeeze(params(i,4,:)), 'Color', colors(i,:), 'LineWidth', 1)
end
xlabel('iteration')
ylabel('\epsilon_\infty')
set(gca,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[3,1,1])

end