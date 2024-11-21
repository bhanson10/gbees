% plot_Lorenz6D.m, https://github.com/bhanson10/gbees-hash/tree/main/examples/Lorenz6D
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

% MC
NM = 1; 
P_DIR = "./results/mc";

count = 1; alpha = 0.05; 
for nm=0:NM-1

    P_DIR_SUB = P_DIR + "/P" + num2str(nm); 
    FILE_LIST = dir(fullfile(P_DIR_SUB, '*.txt'));  % List only .txt files
    num_files = numel(FILE_LIST);

    for i=[0,num_files - 1]
        P_FILE = P_DIR_SUB + "/pdf_" + num2str(i) + ".txt";

        [x_mc, P_mc, n_mc, t_mc(count)] = parse_nongaussian_txt(P_FILE);

        nexttile(1); 
        scatter3(x_mc(:,1), x_mc(:,2), x_mc(:,3), 10, 'filled', 'k', 'HandleVisibility', 'off', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
        nexttile(2);  
        scatter3(x_mc(:,4), x_mc(:,5), x_mc(:,6), 10, 'filled', 'k', 'HandleVisibility', 'off', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
        drawnow;

        count = count + 1;
    end 
end

% %% GBEES
% NM = 1; 
% p.color = "cyan"; p.alpha = [0.3, 0.5, 0.7]; 
% P_DIR = "./results/gbees";
% 
% count = 1;
% for nm=0:NM-1
% 
%     P_DIR_SUB = P_DIR + "/P" + num2str(nm); 
%     FILE_LIST = dir(fullfile(P_DIR_SUB, '*.txt'));  % List only .txt files
%     num_files = numel(FILE_LIST);
% 
%     for i=[0:num_files - 1]
%         P_FILE = P_DIR_SUB + "/pdf_" + num2str(i) + ".txt";
% 
%         [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(P_FILE);
% 
%         xest_gbees{count} = zeros(size(x_gbees(1,:)));
%         for j=1:n_gbees
%             xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
%         end
% 
%         nexttile(1); 
%         plot_nongaussian_surface(x_gbees(:,1:3),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
%         nexttile(2);  
%         plot_nongaussian_surface(x_gbees(:,4:6),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
%         drawnow; 
% 
%         count = count + 1;
%     end
% end

vw0 = [30, 10]; vw = vw0; 
nf = 200; 
dvw = [360 360] ./ nf; 
frames(1) = getframe(gcf); 
for i = 1:nf
    vw = vw + dvw;
    nexttile(1); view(vw(1), vw0(2)); 
    nexttile(2); view(vw(1), vw0(2)); 
    drawnow; 
    pause(0.1); 
    frames(i+1) = getframe(gcf); 
end
create_video(frames,'/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Papers/gbees-gpu/figures/Lorenz96_animation.mp4', 24)

% clear L; clear LH; 
% LH(1) = fill(nan, nan, nan, 'FaceAlpha', 0.7, 'FaceColor', 'cyan', 'EdgeColor', 'none');
% L{1} = "$p_\mathbf{x}(\mathbf{x}', t_{0+})\,\,\,$";
% LH(2) = fill(nan, nan, nan, 'FaceAlpha', 0.7, 'FaceColor', 'magenta', 'EdgeColor', 'none');
% L{2} = "$p_\mathbf{x}(\mathbf{x}', t_{1+})\,\,\,$";
% leg = legend(LH, L, 'Orientation', 'Horizontal', 'FontSize', 18, 'FontName', 'times', 'Interpreter', 'latex');
% leg.Layout.Tile = 'south';
% drawnow; 
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

    f1 = figure(1); clf; hold all; f1.Position = [150 200 1200 475];
    tiledlayout(1, 2, 'TileSpacing','compact');

    nexttile(1); hold all; 
    view(30, 10); lighting phong; light('Position',[-1 0 0]); 
    set(gca, 'FontName' , 'Times','FontSize',12);
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
    
    nexttile(2); hold all; 
    view(30,10); lighting phong; light('Position',[-1 0 0]);
    set(gca, 'FontName' , 'Times','FontSize',12);
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
    P = pdf(:,1); x = pdf(:,2:7); n = size(P, 1);
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%