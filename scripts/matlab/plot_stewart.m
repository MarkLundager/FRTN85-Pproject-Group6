% plot_stewart_log.m
% Real-time serial plot of IMU and REF roll/pitch (Linux)

clear; clc;
 
% --- Porta serial ---
port = '/dev/ttyACM1';     % ajusta conforme necessário
baud = 115200;
s = serialport(port, baud);
configureTerminator(s, "LF");
flush(s);

% --- Inicialização dos vetores ---
t0 = tic;
time = [];
imu_r = [];
imu_p = [];
ref_r = [];
ref_p = [];

figure('Name','Stewart Platform Live Data','NumberTitle','off');
tiledlayout(2,1);

nexttile(1);
hold on;
h1 = plot(nan, nan, 'DisplayName', 'IMU Roll');
h2 = plot(nan, nan, 'DisplayName', 'REF Roll');
yline( 0.5, '--k', 'HandleVisibility','off');
yline(-0.5, '--k', 'HandleVisibility','off');
ylabel('Roll [deg]');
grid on;
legend;
title('Roll');

nexttile(2);
hold on;
h3 = plot(nan, nan, 'DisplayName', 'IMU Pitch');
h4 = plot(nan, nan, 'DisplayName', 'REF Pitch');
yline( 0.5, '--k', 'HandleVisibility','off');
yline(-0.5, '--k', 'HandleVisibility','off');
ylabel('Pitch [deg]');
xlabel('Time [s]');
grid on;
legend;
title('Pitch');

disp('A ler dados... (Ctrl+C para parar)');

while true
    if s.NumBytesAvailable > 0
        line = readline(s);
        tokens = regexp(line, ...
            'IMU r=([-+]?\d*\.?\d+) p=([-+]?\d*\.?\d+).*REF r=([-+]?\d*\.?\d+) p=([-+]?\d*\.?\d+)', ...
            'tokens');
        if ~isempty(tokens)
            vals = str2double(tokens{1});
            t = toc(t0);
            time  = [time,  t];
            imu_r = [imu_r, vals(1)];
            imu_p = [imu_p, vals(2)];
            ref_r = [ref_r, -vals(3)];
            ref_p = [ref_p, -vals(4)];

            set(h1, 'XData', time, 'YData', imu_r);
            set(h2, 'XData', time, 'YData', ref_r);
            set(h3, 'XData', time, 'YData', imu_p);
            set(h4, 'XData', time, 'YData', ref_p);

            drawnow limitrate nocallbacks;
        end
    end
end
