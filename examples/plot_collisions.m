% plot_collisions.m, https://github.com/bhanson10/gbees-hash/tree/main/examples
% Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc; 

CAP = 1024; 
NM = 2; 
C_DIR = "./Lorenz3D/results/c";

figure(1); clf; hold on; axis("equal"); 
xlim([0.5, sqrt(CAP) + 0.5]);
ylim([0.5, sqrt(CAP) + 0.5]);
colormap("abyss");

count = 1;
for nm=0:NM-1

    C_DIR_SUB = C_DIR + "/C" + num2str(nm); 
    FILE_LIST = dir(fullfile(C_DIR_SUB, '*.txt'));  % List only .txt files
    num_files = numel(FILE_LIST);

    for i=[0:num_files - 1]
        C_FILE = C_DIR_SUB + "/c_" + num2str(i) + ".txt";

        c = parse_collisions_files(C_FILE);
        clims = [min(c), max(c)];
        n = sqrt(size(c, 1)); 
        c = reshape(c, [n, n]); 
        imagesc(c, clims); 
        colorbar; 
        title("Step: " + num2str(nm) + "-" + num2str(i), 'FontSize', 16, 'FontName', 'times'); 
        drawnow; 

        frames(count) = getframe(gcf); 
        count = count + 1; 
    end
end

% create_video(frames, "Lorenz3D_collisions.mp4", 20);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = parse_collisions_files(FILE_NAME)
    c = []; 
    f = fopen(FILE_NAME, 'r');
    while ~feof(f)
        line = split(fgetl(f)); 
        c(end+1, 1) = str2double(line{1});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%