% plot_duffing.m, https://github.com/bhanson10/gbees/tree/main/examples/duffing
% Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc; format shortG;

%% parameters
load("colors.mat"); 
mu_s = [0; 1]; S_s = [0.01 0; 0 0.01]; N = 2; 

% trajectory
tspan = [0 100]; % propagate from t=0 to t=100
[~, Xs] = ode45(@(t,x) duffing(t,x), tspan, mu_s);

% mc
Nmc = 5000; 
X0 = mvnrnd(mu_s, S_s, Nmc); 

% CUT
[Z0, W] = cut6(mu_s, S_s); n = size(Z0, 1); 

% plotting
figure("Position", [348,206,784,604]); 
tiledlayout(1, 1, "TileSpacing", "compact"); nexttile(1); hold on; box on;
set(gca, 'FontName', 'times', 'FontSize', 18, "LineWidth", 2);  
scatter(X0(:, 1), X0(:, 2), 50, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerFaceAlpha', 0.25);
scatter(Z0(:, 1), Z0(:, 2), 100, 'd', 'filled', 'MarkerFaceColor', hanred);

% Principle of maximum entropy
M = 6; % up to 4th order
[X, P] = pme_cut2pdf(Z0, W, 'M', M); 

% GBEES
P_FILE = "./results/c/P0/pdf_0.txt";
[x_gbees, P_gbees, ~, ~] = parse_nongaussian_txt(P_FILE);

p.lw = 3; p.ls = "-"; p.color = hangreen; 
plot_nongaussian_surface(X, 'P', P, 'p', p);
plot(Xs(:, 1), Xs(:, 2), "g-", "LineWidth", 1);
p.ls = "-."; p.color = hanblue; 
plot_nongaussian_surface(x_gbees, 'P', P_gbees, 'p', p); 
drawnow; 

for tend = [1, 2.5, 5]
    tspan = [0 tend]; 
    [~, Xt] = ode45(@(t,x) duffing(t,x), tspan, mu_s);
    plot(Xt(:, 1), Xt(:, 2), "k-", "LineWidth", 2);

    Xf = zeros(Nmc, 2); % store propagated states
    for i = 1:Nmc
        [~, Xi] = ode45(@(t,x) duffing(t,x), tspan, X0(i,:)');
        Xf(i, :) = Xi(end, :); 
    end
    scatter(Xf(:, 1), Xf(:, 2), 50, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerFaceAlpha', 0.25);
    
    Zf = zeros(n, 2); % store propagated states
    for i = 1:n
        [~, Zi] = ode45(@(t,x) duffing(t,x), tspan, Z0(i,:)');
        Zf(i, :) = Zi(end, :); 
    end
    scatter(Zf(:, 1), Zf(:, 2), 100, 'd', 'filled', 'MarkerFaceColor', hanred);

    % Principle of Maximum Entropy
    mu = reconstruct_cut(Zf, W, 1, 4); % 1) first central moment (mean)   
    S  = reconstruct_cut(Zf, W, 2, 4); % 2) second central moment (covariance)
    [X, P] = pme_cut2pdf(Zf, W, 'M', M); 

    % GBEES
    file_num = (tend / 0.5); 
    P_FILE = "./results/<langugage>/P0/pdf_" + num2str(file_num) + ".txt";
    [x_gbees, P_gbees, ~, ~] = parse_nongaussian_txt(P_FILE);

    p.ls = "-"; p.color = hangreen; 
    plot_nongaussian_surface(X, 'P', P, 'p', p);
    p.ls = "-."; p.color = hanblue; 
    plot_nongaussian_surface(x_gbees, 'P', P_gbees, 'p', p);
    drawnow; 
end

% legend
LH(1) = scatter(nan, nan, 100, 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
L{1} = "MC $(n=$ " + num2str(Nmc) + "$)$";
LH(2) = scatter(nan, nan, 100, 'filled', 'd', 'filled', 'MarkerFaceColor', hanred, 'MarkerEdgeColor', 'none');
L{2} = "CUT6 $(N=$ " + num2str(numel(W)) + "$)$";
LH(3) = plot(nan, nan,  "Color", hangreen, "LineWidth", 2);
L{3} = "PME ($M=6$)";
LH(4) = plot(nan, nan,  "Color", hanblue, "LineWidth", 2);
L{4} = "GBEES";
leg = legend(LH, L, "FontSize", 18, "FontName", "times", "Interpreter", "latex",...
             "LineWidth", 1, "Orientation", "horizontal");
leg.Layout.Tile = 'south'; 
xlim([-1.75, 1.75]); ylim([-1.25, 1.5]); 

text(-0.13671875,0.5369211514393, "$t=0$", "FontSize", 20, "Interpreter", "latex");
text(0.8818359375,0.764080100125156, "$t=1$", "FontSize", 20, "Interpreter", "latex");
text(0.765625,-0.619524405506883, "$t=2.5$", "FontSize", 20, "Interpreter", "latex");
text(-0.8203125,-0.998122653316646, "$t=5$", "FontSize", 20, "Interpreter", "latex");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = duffing(t, x)
    a = 0.35; b = 0.3; w = 1; 
    dx = zeros(2, 1); 
    dx(1) = x(2); 
    dx(2) = x(1) - x(1)^3 - a * x(2) + b * cos(w * t); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename)
    fileID = fopen(filename, 'r');
    data = textscan(fileID, '%s', 'Delimiter', '\n'); data = data{1};
    t = str2num(data{1}); data = data(2:end); 
    pdf = cellfun(@(x) str2num(x), data, 'UniformOutput', false); pdf = cell2mat(pdf); 
    P = pdf(:,1); x = pdf(:,2:end); n = size(P, 1);
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%