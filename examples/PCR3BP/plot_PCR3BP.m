% plot_PCR3BP.m, https://github.com/bhanson10/gbees/tree/main/examples/PCR3BP
% Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc;

%% initializing system properties
prop.mu = 1.901109735892602e-07; prop.LU = 238529; prop.TU = 18913; prop.sec_r = 2.521000000000000e+02;
prop.U = [prop.LU, prop.LU, prop.LU/prop.TU, prop.LU/prop.TU]; 
ic.T = 3.727168019157753;
ic.state = [1.001471995170839 -0.000017518099335 0.000071987832396 0.013633926328993];

%% initializing figure
f1 = figure(1); clf; hold all; f1.Position = [250 100 1000 700];
tiledlayout(1, 2, 'TileSpacing','compact');

nexttile(1); hold all; 
axis("equal")
set(gca, 'FontName', 'Times', 'FontSize', 14);
xlabel("synodic $x$ (km)", "Interpreter","latex")
ylabel("synodic $y$ (km)", "Interpreter","latex")
europa = nsidedpoly(1000, 'Center', [(1-prop.mu)*prop.LU, 0], 'Radius', prop.sec_r);
plot(europa, 'FaceColor', 'm', 'EdgeColor','none');

nexttile(2); hold all; 
axis("equal")
set(gca, 'FontName', 'Times', 'FontSize', 14);
xlabel("synodic $\dot{x}$ (km/s)", "Interpreter","latex")
ylabel("synodic $\dot{y}$ (km/s)", "Interpreter","latex")
 
%% truth
options = odeset('MaxStep', 1E-3, 'InitialStep', 1E-3, 'RelTol', 1e-6);
[t, x] = ode87(@(t, x) PCR3BP(t, x, prop), [0,ic.T], ic.state, options);
x(:,1:4) = x(:,1:4).*prop.U; 
t = t.*prop.TU; 

%% nominal trajectory
nexttile(1); 
plot(x(:,1),x(:,2),'k-','LineWidth', 2,'HandleVisibility','off');
drawnow;
nexttile(2); 
plot(x(:,3),x(:,4),'k-','LineWidth', 2,'HandleVisibility','off');
drawnow;

%% GBEES
NM = 4; 

C = [238, 102, 119;  % Red
     68,  119, 170;  % Blue
     255, 140, 0;    % Orange
     34,  139, 34;]; % Green

C  = C/255;

p.alpha = [0.3, 0.4, 0.7]; 
P_DIR = "<path_to_pdf>";

count = 1;
for nm=0:NM-1
    p.color = C(nm + 1, :);
    P_DIR_SUB = P_DIR + "/P" + num2str(nm); 
    FILE_LIST = dir(fullfile(P_DIR_SUB, '*.txt'));  % List only .txt files
    num_files = numel(FILE_LIST);
    
    if mod(nm, 2) == 0 % alternate between blue and red isocurves
        set = [0, 2, num_files-1];
    else
        set = [0, 5, num_files-1];
    end

    for i=set
        P_FILE = P_DIR_SUB + "/pdf_" + num2str(i) + ".txt";

        [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(P_FILE, prop);
    
        xest_gbees{count} = zeros(size(x_gbees(1,:)));
        for j=1:n_gbees
            xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
        end

        nexttile(1); 
        plot_nongaussian_surface(x_gbees(:,1:2),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
        nexttile(2); 
        plot_nongaussian_surface(x_gbees(:,3:4),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)], p);
        drawnow; 
        
        count = count + 1;
    end
end

P_FILE = P_DIR + "/P0/pdf_0.txt";
p.color = C(1, :);
[x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(P_FILE, prop);
nexttile(1); 
plot_nongaussian_surface(x_gbees(:,1:2),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
nexttile(2); 
plot_nongaussian_surface(x_gbees(:,3:4),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)], p);
drawnow; 

clear L; clear LH; 
LH(1) = plot(NaN,NaN,'o','MarkerSize', 10, 'MarkerFaceColor',[1,167/255,254/255],'MarkerEdgeColor','none');
L{1} = "Enceladus\,\,\,";
LH(2) = plot(NaN,NaN,'k-', 'LineWidth',1);
L{2} = "Nominal\,\,\,";
LH(3) = fill(nan, nan, nan, 'FaceAlpha', 0.7, 'FaceColor', C(1,:), 'EdgeColor', 'none');
L{3} = "$p_\mathbf{x}(\mathbf{x}', t_{0+})\,\,\,$";
LH(4) = fill(nan, nan, nan, 'FaceAlpha', 0.7, 'FaceColor', C(2,:), 'EdgeColor', 'none');
L{4} = "$p_\mathbf{x}(\mathbf{x}', t_{1+})\,\,\,$";
LH(5) = fill(nan, nan, nan, 'FaceAlpha', 0.7, 'FaceColor', C(3,:), 'EdgeColor', 'none');
L{5} = "$p_\mathbf{x}(\mathbf{x}', t_{2+})\,\,\,$";
LH(6) = fill(nan, nan, nan, 'FaceAlpha', 0.7, 'FaceColor', C(4,:)', 'EdgeColor', 'none');
L{6} = "$p_\mathbf{x}(\mathbf{x}', t_{3+})\,\,\,$";
leg = legend(LH, L, 'Orientation', 'Horizontal', 'FontSize', 18, 'FontName', 'times', 'Interpreter', 'latex');
leg.Layout.Tile = 'south';
drawnow; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = PCR3BP(t, x, prop)
    x1 = [x(3); x(4); 2*x(4)+x(1)-(prop.mu*(x(1)-1+prop.mu)/(((x(1)-1+prop.mu)^2+x(2)^2)^(1.5)))-((1-prop.mu)*(x(1)+prop.mu)/(((x(1)+prop.mu)^2+x(2)^2)^(1.5))); -2*x(3)+x(2)-(prop.mu*x(2)/(((x(1)-1+prop.mu)^2+x(2)^2)^(1.5)))-((1-prop.mu)*x(2)/(((x(1)+prop.mu)^2+x(2)^2)^(1.5)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename, prop)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID));
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count,1) = str2double(line{1});
        x(count, :) = [str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5})].*prop.U';
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%