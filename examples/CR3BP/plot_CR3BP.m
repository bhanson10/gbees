% plot_CR3BP.m, https://github.com/bhanson10/gbees/tree/main/examples/CR3BP
% Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc;

%% Initial Condition
prop.mu = 2.528017528540000E-5; prop.LU = 668519; prop.TU = 48562; prop.sec_r =  1560.8;    
prop.U = [prop.LU, prop.LU, prop.LU, prop.LU/prop.TU, prop.LU/prop.TU, prop.LU/prop.TU]; 
prop.T = 2.6513344042156235E+0; prop.r = 1560.8;
ic = [1.0169962521469986; -1.0697953862098174E-20; -5.1360139986508944E-34; -1.3935178510137174E-14; 1.2575912545456968E-2; -3.1572204202682494E-33];

%% initializing figure
f1 = figure(1); clf; hold on; f1.Position = [50 100 1400 700];
tiledlayout(1, 2, 'TileSpacing','compact');

nexttile(1); hold on; 
axis("equal")
set(gca, 'FontName', 'Times', 'FontSize', 14);
xlabel("synodic $x$ (km)", "Interpreter","latex")
ylabel("synodic $y$ (km)", "Interpreter","latex")
zlabel("synodic $z$ (km)", "Interpreter","latex")
[X, Y, Z] = sphere(100);  
X = X.*prop.r + (1-prop.mu)*prop.LU; Y = Y.*prop.r; Z = Z.*prop.r;
surf(X, Y, Z, 'FaceColor', 'g', 'EdgeColor', 'none'); alpha(0.3); 
drawnow;

nexttile(2); hold on; 
axis("equal")
set(gca, 'FontName', 'Times', 'FontSize', 14);
xlabel("synodic $\dot{x}$ (km/s)", "Interpreter","latex")
ylabel("synodic $\dot{y}$ (km/s)", "Interpreter","latex")
zlabel("synodic $\dot{z}$ (km/s)", "Interpreter","latex")

%% truth
options = odeset('MaxStep', 1E-3, 'InitialStep', 1E-3, 'RelTol', 1e-6);
[t, x] = ode87(@(t, x) CR3BP(t, x, prop), [0, prop.T], ic, options);
x(:,1:6) = x(:,1:6).*prop.U; 
t = t.*prop.TU; 

%% nominal trajectory
nexttile(1); 
plot3(x((t <= 12*3600),1),x((t <= 12*3600),2),x((t <= 12*3600),3),'r-','LineWidth',2,'DisplayName','Nominal');
plot3(x(:,1),x(:,2),x(:,3),'r--','LineWidth',1,'DisplayName','Nominal');
drawnow;
nexttile(2); 
plot3(x((t <= 12*3600),4),x((t <= 12*3600),5),x((t <= 12*3600),6),'r-','LineWidth',2,'DisplayName','Nominal');
plot3(x(:,4),x(:,5),x(:,6),'r--','LineWidth', 1, 'DisplayName','Nominal');
drawnow;

%% GBEES
binary = 1;
if binary, file_str = ".bin"; else, file_str = ".txt"; end 
NM = 1; 
p.color = 'b'; p.type = "grid"; 
P_DIR = "./results/<language>";

count = 1;
for nm=0:NM-1

    P_DIR_SUB = P_DIR + "/P" + num2str(nm); 
    FILE_LIST = dir(fullfile(P_DIR_SUB, '*' + file_str));  % List only file_str files
    num_files = numel(FILE_LIST);

    for i=[0, num_files - 1]
        P_FILE = P_DIR_SUB + "/pdf_" + num2str(i) + file_str;

        [x_gbees, P_gbees, ~, ~] = parse_nongaussian_txt(P_FILE, binary, prop);

        nexttile(1); 
        plot_nongaussian_surface(x_gbees(:,1:3), 'P', P_gbees, 'p', p);
        drawnow; 
        nexttile(2); 
        plot_nongaussian_surface(x_gbees(:,4:6), 'P', P_gbees, 'p', p);
        drawnow; 
        count = count + 1;
    end
end

clear L; clear LH; 
LH(1) = plot(NaN,NaN,"o","MarkerSize", 10, "MarkerFaceColor","g","MarkerEdgeColor","k");
L{1} = "Europa\,\,\,";
LH(2) = plot(NaN,NaN,"k-", "LineWidth",1);
L{2} = "Nominal\,\,\,";
LH(3) = fill(nan, nan, nan, "FaceAlpha", 0.7, "FaceColor", 'b', "EdgeColor", "none");
L{3} = "$p(\mathbf{x},t = [0,12 \text{ hours}])\,\,\,$";
leg = legend(gca, LH, L, "Orientation", "Horizontal", "Position", [0.374216319492885, 0.01, 0.25156736101423, 0.075739644688262], "FontSize", 18, "FontName", "times", "Interpreter", "latex");
leg.Layout.Tile = "south";
drawnow; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xk = CR3BP(t, x, prop)
    xk = [x(4); x(5); x(6); 2*x(5)+x(1)-(prop.mu*(x(1)-1+prop.mu)/(((x(1)-1+prop.mu)^2+x(2)^2+x(3)^2)^(1.5)))-((1-prop.mu)*(x(1)+prop.mu)/(((x(1)+prop.mu)^2+x(2)^2+x(3)^2)^(1.5))); -2*x(4)+x(2)-(prop.mu*x(2)/(((x(1)-1+prop.mu)^2+x(2)^2+x(3)^2)^(1.5)))-((1-prop.mu)*x(2)/(((x(1)+prop.mu)^2+x(2)^2+x(3)^2)^(1.5))); -(prop.mu*x(3)/(((x(1)-1+prop.mu)^2+x(2)^2+x(3)^2)^(1.5)))-((1-prop.mu)*x(3)/(((x(1)+prop.mu)^2+x(2)^2+x(3)^2)^(1.5)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename, prop, binary)
    if binary
        fid = fopen(filename, 'rb');
        assert(fid ~= -1, 'Could not open file.');
        t = fread(fid, 1, 'double');
        n = fread(fid, 1, 'uint32');
        raw = fread(fid, [7, n], 'double');
        P = raw(1,:)';
        x = raw(2:end,:)'.*prop.U;
        fclose(fid);
    else
        fileID = fopen(filename, 'r');
        data = textscan(fileID, '%s', 'Delimiter', '\n'); data = data{1};
        t = str2num(data{1}); data = data(2:end); 
        pdf = cellfun(@(x) str2num(x), data, 'UniformOutput', false); pdf = cell2mat(pdf); 
        P = pdf(:,1); x = pdf(:,2:end).*prop.U; n = size(P, 1);
        fclose(fileID);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%