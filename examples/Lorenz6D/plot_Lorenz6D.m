% plot_Lorenz6D.m, https://github.com/bhanson10/gbees/tree/main/examples/Lorenz6D
% Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

close all; clc; clear all; 

%% initializing system properties
prop.T = 1.5; prop.F = 4;
ic = [prop.F + 0.5; prop.F; prop.F; prop.F; prop.F; prop.F];

%% initializing figure
initialize_figures(); 

%% truth
options = odeset('MaxStep', 1E-3, 'InitialStep', 1E-3, 'RelTol', 1e-6);
[~, xf] = ode87(@(t, xf) Lorenz6D(t, xf, prop), [0 50], ic, options);
[~, x] = ode87(@(t, x) Lorenz6D(t, x, prop), [0 prop.T], ic, options);

nexttile(1); 
plot3(xf(:,1), xf(:,2), xf(:,3), 'g-','linewidth', 1, 'HandleVisibility', 'off');
plot3(x(:,1), x(:,2), x(:,3), 'k-','linewidth', 2, 'HandleVisibility', 'off');    
drawnow;

nexttile(2); 
plot3(xf(:,4),xf(:,5),xf(:,6),'g-','linewidth', 1,'DisplayName','Nominal');
plot3(x(:,4),x(:,5),x(:,6),'k-','linewidth',2,'DisplayName','Nominal'); 
drawnow;

nexttile(3); 
plot3(xf(:,1), xf(:,2), xf(:,3), 'g-','linewidth', 1, 'HandleVisibility', 'off');
plot3(x(:,1), x(:,2), x(:,3), 'k-','linewidth', 2, 'HandleVisibility', 'off');    
drawnow;

nexttile(4); 
plot3(xf(:,4),xf(:,5),xf(:,6),'g-','linewidth', 1,'DisplayName','Nominal');
plot3(x(:,4),x(:,5),x(:,6),'k-','linewidth',2,'DisplayName','Nominal'); 
drawnow;

% MC
NM = 1; 
P_DIR = "./results/mc";

count = 1; alpha = 0.05; 
for nm=0:NM-1

    FILE_LIST = dir(fullfile(P_DIR, '*.txt'));  % List only .txt files
    num_files = numel(FILE_LIST);

    for i=[0,num_files - 1]
        P_FILE = P_DIR + "/mc_" + num2str(i) + ".txt";

        [x_mc, P_mc, n_mc, t_mc(count)] = parse_nongaussian_txt(P_FILE);

        nexttile(1); 
        scatter3(x_mc(:,1), x_mc(:,2), x_mc(:,3), 10, 'filled', 'k', 'HandleVisibility', 'off', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
        nexttile(2);  
        scatter3(x_mc(:,4), x_mc(:,5), x_mc(:,6), 10, 'filled', 'k', 'HandleVisibility', 'off', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
        drawnow;

        count = count + 1;
    end 
end

%% GBEES
NM = 1; 
p.color = "cyan"; p.type = "grid"; 
P_DIR = "./results/gbees/<language>";

count = 1;
for nm=0:NM-1

    P_DIR_SUB = P_DIR + "/P" + num2str(nm); 
    FILE_LIST = dir(fullfile(P_DIR_SUB, '*.txt'));  % List only .txt files
    num_files = numel(FILE_LIST);

    for i=[0:num_files - 1]
        P_FILE = P_DIR_SUB + "/pdf_" + num2str(i) + ".txt";

        [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(P_FILE);

        [x_list, ~, x_idx] = unique(x_gbees(:,1:3), 'rows', 'stable');
        Px_full = accumarray(x_idx, P_gbees, [], @sum);

        [v_list, ~, v_idx] = unique(x_gbees(:,4:6), 'rows', 'stable');
        Pv_full = accumarray(v_idx, P_gbees, [], @sum);

        nexttile(3); 
        plot_nongaussian_surface(x_list, Px_full, 'p', p);
        nexttile(4);  
        plot_nongaussian_surface(v_list, Pv_full, 'p', p);
        axis normal
        drawnow; 

        count = count + 1;
    end
end

% create video
% vw0 = [30, 10]; vw = vw0; 
% nf = 200; 
% dvw = [360 360] ./ nf; 
% frames(1) = getframe(gcf); 
% for i = 1:nf
%     vw = vw + dvw;
%     nexttile(1); view(vw(1), vw0(2)); 
%     nexttile(2); view(vw(1), vw0(2)); 
%     drawnow; 
%     pause(0.1); 
%     frames(i+1) = getframe(gcf); 
% end
% create_video(frames,'Lorenz96_animation.mp4', 24)

clear L; clear LH; 
LH(1) = scatter(nan, nan, 10, 'k', 'filled');
L{1} = "MC {      }";
LH(2) = fill(nan, nan, nan, 'FaceAlpha', 0.7, 'FaceColor', 'cyan', 'EdgeColor', 'none');
L{2} = "$p(\mathbf{x}, t = [0,1.3])\,\,\,$";
leg = legend(LH, L, 'Orientation', 'Horizontal', 'FontSize', 18, 'FontName', 'times', 'Interpreter', 'latex');
leg.Layout.Tile = 'south';
drawnow; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              FUNCTIONS                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=Lorenz6D(t, x, prop)                          
    f = [(x(2) - x(5)) * x(6) - x(1) + prop.F;
         (x(3) - x(6)) * x(1) - x(2) + prop.F;
         (x(4) - x(1)) * x(2) - x(3) + prop.F;
         (x(5) - x(2)) * x(3) - x(4) + prop.F;
         (x(6) - x(3)) * x(4) - x(5) + prop.F;
         (x(1) - x(4)) * x(5) - x(6) + prop.F];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures()

    f1 = figure(1); clf; hold on; f1.Position = [100 50 1200 850];
    tiledlayout(2, 2, 'TileSpacing','compact');

    nexttile(1); hold on; axis normal
    view(30, 10); lighting phong; light('Position',[1 -1 1]); 
    set(gca, 'FontName' , 'Times','FontSize',14);
    xlabel("$x_1$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel("$x_2$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    zlabel("$x_3$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    set(get(gca,'ZLabel'), 'Rotation', 0);
    xlim([-6 8])
    xticks([-6 -4 -2 0 2 4 6 8])
    xticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    ylim([-6 8])
    yticks([-6 -4 -2 0 2 4 6 8])
    yticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    zlim([-4 8])
    zticks([-4 -2 0 2 4 6 8])
    zticklabels({'-4','-2','0','2','4', '6', '8'})
    set(gca, 'FontName' , 'Times');
    text(-13.5, 0, 1, '(a)', 'FontSize', 30, 'FontName', 'times', 'Interpreter', 'latex');
    
    nexttile(2); hold on; axis normal
    view(30,10); lighting phong; light('Position',[1 -1 1]);
    set(gca, 'FontName' , 'Times','FontSize',14);
    xlabel("$x_4$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel("$x_5$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    zlabel("$x_6$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    set(get(gca,'ZLabel'), 'Rotation', 0);
    xlim([-6 8])
    xticks([-6 -4 -2 0 2 4 6 8])
    xticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    ylim([-6 8])
    yticks([-6 -4 -2 0 2 4 6 8])
    yticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    zlim([-6 8])
    zticks([-6 -4 -2 0 2 4 6 8])
    zticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    set(gca, 'FontName' , 'Times');

    nexttile(3); hold on; axis normal
    view(30, 10); lighting phong; light('Position',[1 -1 1]); 
    set(gca, 'FontName' , 'Times','FontSize',14);
    xlabel("$x_1$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel("$x_2$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    zlabel("$x_3$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    set(get(gca,'ZLabel'), 'Rotation', 0);
    xlim([-6 8])
    xticks([-6 -4 -2 0 2 4 6 8])
    xticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    ylim([-6 8])
    yticks([-6 -4 -2 0 2 4 6 8])
    yticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    zlim([-4 8])
    zticks([-4 -2 0 2 4 6 8])
    zticklabels({'-4','-2','0','2','4', '6', '8'})
    set(gca, 'FontName' , 'Times');
    text(-13.5, 0, 1, '(b)', 'FontSize', 30, 'FontName', 'times', 'Interpreter', 'latex');
    
    nexttile(4); hold on; axis normal 
    view(30,10); lighting phong; light('Position',[1 -1 1]);
    set(gca, 'FontName' , 'Times','FontSize',14);
    xlabel("$x_4$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel("$x_5$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    zlabel("$x_6$", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
    set(get(gca,'ZLabel'), 'Rotation', 0);
    xlim([-6 8])
    xticks([-6 -4 -2 0 2 4 6 8])
    xticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    ylim([-6 8])
    yticks([-6 -4 -2 0 2 4 6 8])
    yticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    zlim([-6 8])
    zticks([-6 -4 -2 0 2 4 6 8])
    zticklabels({'-6', '-4','-2','0','2','4', '6', '8'})
    set(gca, 'FontName' , 'Times');
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