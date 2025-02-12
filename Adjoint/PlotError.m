function [] = PlotError()

% Filenames
filename = {'phasesep.mat', 'phasesep_T0.05.mat', 'phasesep_T0.10.mat', 'phasesep_T0.20.mat', 'phasesep_T0.50.mat'};
ldg_labels = {'0.00', '0.05', '0.10', '0.20', '0.50'};
lgd_title = {'T'};

% Setup figures
close all
h_error = figure; hold on
colors = parula(length(filename));

for i = 1:length(filename)

    % Load data
    load(filename{i})
    
    M = size(E,1);
    plot(1:M, E, '-', 'Color', colors(i,:), 'LineWidth', 1.5)

end

% Format error plot
figure(h_error)
xlabel('iteration')
ylabel('error')
set(gca, 'XScale', 'log', 'YScale', 'log')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])

% Add legend
ldg = legend(ldg_labels);
title(ldg, lgd_title)

end