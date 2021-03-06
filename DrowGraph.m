clear;
close all;
clc;
% reading file
fileID = fopen('LQR_Kalman.txt', 'r');
formatSpec = '%f';
size = [13 251];
data = fscanf(fileID, formatSpec, size);
% graph
fig1 = figure(1);
set(groot, 'DefaultTextInterpreter', 'Latex');
set(groot, 'DefaultLegendInterpreter', 'Latex');
for i = 1 : 4
    subplot(2, 2, i);
    plot(data(1, :), data(i + 1, :), 'b'); hold on;
    plot(data(1, :), data(i + 9, :), '--r'); hold on;
    xlim([0 25]);
    switch i
        case 1
            xlabel('time[sec]');
            ylabel('$x$');
            legend('$x$');
        case 2
            xlabel('time[sec]');
            ylabel('$\theta$');
            legend('$\theta$');
        case 3
            xlabel('time[sec]');
            ylabel('$\dot{x}$');
            legend('$\dot{x}$');
        case 4
            xlabel('time[sec]');
            ylabel('$\dot{\theta}$');
            legend('$\dot{\theta}$');
        otherwise
            disp('Error');
    end
    grid on;
end


fig2 = figure(2);
set(groot, 'DefaultTextInterpreter', 'Latex');
set(groot, 'DefaultLegendInterpreter', 'Latex');
for i = 1 : 4
    subplot(2, 2, i);
    plot(data(1, :), data(i + 1, :), 'b'); hold on;
    plot(data(1, :), data(i + 5, :), '--r'); hold on;
    xlim([0 25]);
    switch i
        case 1
            xlabel('time[sec]');
            ylabel('$x$');
            legend('$x$');
        case 2
            xlabel('time[sec]');
            ylabel('$\theta$');
            legend('$\theta$');
        case 3
            xlabel('time[sec]');
            ylabel('$\dot{x}$');
            legend('$\dot{x}$');
        case 4
            xlabel('time[sec]');
            ylabel('$\dot{\theta}$');
            legend('$\dot{\theta}$');
        otherwise
            disp('Error');
    end
    grid on;
end

% Auto save graph
saveas(fig1, 'LQR_result.png');
saveas(fig2, 'Estimate.png');