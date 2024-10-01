% plot_collisions.m, https://github.com/bhanson10/gbees-hash/tree/main/examples
% Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc; 

CAP = 2^10; 
NM = 2; 
SYS = "Lorenz3D";
LANG = "c"; 
HASH = "XNOR";

n = sqrt(CAP); 
figure(1); clf; hold on; axis("square"); 
set(gca, "FontSize", 14, "FontName", "times");
xlim([0.5, n + 0.5]);
ylim([0.5, n + 0.5]);
colormap("abyss");

count = 1;
for nm=0:NM-1

    C_DIR = "./" + SYS + "/results/" + LANG + "/C" + num2str(nm); 
    FILE_LIST = dir(fullfile(C_DIR, '*.txt'));  % List only .txt files
    num_files = numel(FILE_LIST);

    for i=[0:num_files - 1]
        C_FILE = C_DIR + "/c_" + num2str(i) + ".txt";

        c = parse_collisions_files(C_FILE);
        clims = [min(c), max(c)];
        c = reshape(c, [n, n]); 
        imagesc(c, clims); 
        cb = colorbar;
        ylabel(cb, "# of collisions at entry", 'FontSize', 18, 'FontName', 'times');
        title(HASH + ",   Capacity: " + CAP + ",   Step: " + num2str(nm) + "-" + num2str(i), 'FontSize', 16, 'FontName', 'times'); 
        drawnow; 

        frames(count) = getframe(gcf); 
        count = count + 1; 
        return
    end
end

% create_video(frames, "../../figures/" + SYS + "/" + SYS + "_" + HASH + "_collisions.mp4", 20);
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