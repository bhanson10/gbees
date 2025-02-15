% plot_POE2BP.m, https://github.com/bhanson10/gbees/tree/main/examples/POE2BP
% Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

close all; clc; clear all; 

%% initializing figure
initialize_figures(); 

%% MC
x0 = [4.6679; 0]; P0 = [4.4723E-5 0.0000000; 0.0000000 3.0461E-8]; N = 1e3; 
X0 = mvnrnd(x0, P0, N); t = linspace(0, 32.242, 5); 
x_mc_1 = zeros(N, size(t,2), 2);
for i = 1:N
    [~, X] = ode45(@(t, x) R2BP(t, x), t, X0(i, :)');
    
    for j = 1:5
        x_mc_1(i, j, :) = X(j, :);
    end
end

P0 = [0.06243 0.0000000; 0.0000000 3.0461E-8];  
X0 = mvnrnd(x0, P0, N); t = linspace(0, 32.242, 5); 
x_mc_2 = zeros(N, size(t,2), 2);
for i = 1:N
    [~, X] = ode45(@(t, x) R2BP(t, x), t, X0(i, :)');
    
    for j = 1:5
        x_mc_2(i, j, :) = X(j, :);
    end
end

count = 1; 
p1.type = "line"; p1.color = "r"; p1.display = 0; p1.name = "$3\sigma$";
p2.color = "jet"; p2.display = 0; p2.name = "GBEES";
for i = 1:size(x_mc_1, 2)
    nexttile(count); title(num2str(5 * (i-1)) + " orbits", 'FontSize', 18, 'FontName', 'Times'); 
    if i ~= 1, xlim([-2, 2]); ylim([-2e-4, 1e-4]); else, xlim([-0.03, 0.03]); ylim([-8e-4, 8e-4]); end
    if i == size(x_mc_1, 2)
        leg = legend("Color", "w", "Orientation", "Horizontal", "FontSize", 24, "FontName", "times", "Interpreter", "latex");
        leg.Layout.Tile = "south";
        p1.display = 1; 
    end

    % MC
    x_mc = [x_mc_1(:, i, 1) x_mc_1(:, i, 2)];
    xi = mean(x_mc); 
    Pi = cov(x_mc);
    [Vi, Di] = eig(Pi); R = Vi * [0 1; -1 0];
    x_mc = (x_mc - xi) * R;
    
    % GBEES
    P_FILE = "./results/<language>/P0/pdf_" + num2str(i - 1) + ".txt";
    [x_gbees, P_gbees, n_gbees, ~] = parse_nongaussian_txt(P_FILE);
    
    Nx = size(unique(x_gbees(:,1)),1); Ny = size(unique(x_gbees(:,2)),1);
    x_gbees = (x_gbees - xi) * R; 
   
    % plotting  
    scatter(x_gbees(:, 1), x_gbees(:, 2), 10, P_gbees, "filled", "DisplayName", "GBEES");  
    scatter(x_mc(:, 1), x_mc(:, 2), 20, "r", "Marker", "x", "DisplayName", "MC {     }");   
    plot_gaussian_ellipsoid([0; 0], R * Pi * R', 3, p1);
    count = count + 1; 
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              FUNCTIONS                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures()

    cmap = jet; bg = cmap(1, :); 

    f1 = figure(1); clf; hold on; f1.Position = [30,232,1461,372];
    tiledlayout(1, 5, 'TileSpacing','compact');
    
    count = 1;
    for i = 1:1
        for j = 1:5
            nexttile(count); hold on; set(gca, 'Color', bg, 'FontName' , 'Times', 'FontSize', 14); colormap("jet");
            if (i == 1)&&(j==1)
                xlabel("$\delta L$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex'); 
                ylabel("$\delta l$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex'); 
            end

            count = count + 1; 
        end
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = R2BP(t, x)    
    mu = 19.910350621818949;
    xdot = [0; mu^2 / x(1)^3];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename)
    fileID = fopen(filename, 'r');
    data = textscan(fileID, '%s', 'Delimiter', '\n'); data = data{1};
    t = str2num(data{1}); data = data(2:end); 
    pdf = cellfun(@(x) str2num(x), data, 'UniformOutput', false); pdf = cell2mat(pdf); 
    P = pdf(:,1); x = pdf(:,2:3); n = size(P, 1);
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%