% compare_times.m, https://github.com/bhanson10/gbees/tree/main/examples
% Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc; 

mkrs = {"-o", "-square", "->", "-pentagram"}; 
SYS = "<system>";
DIR = {"c", "python"}; 
TU = 1; 
names = {'c','python'};
time_unit = '(TU)';

f = {}; 
for i=1:length(names)
    f{end+1} = fopen("./" + SYS + "/results/" + DIR{i} + "/runtime.txt", 'r');
end

[st, pt, s] = parse_time_files(f);

f1 = figure(1); clf; f1.Position = [100 350 1300 475];
tiles = tiledlayout(1,2);

nexttile(1); ax = gca; 
set(ax, 'FontName' , 'Times','FontSize',14);
xlabel("Simulation time " + time_unit, 'FontSize', 18, 'FontName', 'Times');
yyaxis left; hold on; 
ax.YAxis(1).Color = [0 0 1];
for i=1:length(f)
    sti = st{i}; pti = pt{i}; 
    plot(sti.*TU, pti, mkrs{i}, 'Color', 'blue', 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName', names{i});
end
ylabel("Program time (s)", 'FontSize', 18, 'FontName', 'Times');

yyaxis right; hold on; 
ax.YAxis(2).Color = [1 0 0];
for i=1:length(f)
    sti = st{i}; si = s{i}; 
    plot(sti.*TU, si, mkrs{i}, 'Color', 'red', 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName', names{i});
end
ylabel("Cell count", 'FontSize', 18, 'FontName', 'Times');
xlim([0, sti(end).*TU]); 

for i = 1:length(f)
    LH1(i) = plot(NaN, NaN, mkrs{i}, 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 10);
    L1{i} = names{i};
end
leg1 = legend(LH1, L1, 'Location', 'North', 'Orientation', 'Horizontal', 'FontSize', 14, 'FontName', 'times', 'Interpreter', 'latex');
leg1.Layout.Tile = "south";
drawnow; 

nexttile(2); ax = gca; 
set(ax, 'FontName' , 'Times','FontSize',14);
xlabel("Simulation time " + time_unit, 'FontSize', 18, 'FontName', 'Times');

yyaxis left; hold on; 
st1 = st{1}; pt1 = pt{1}; name1 = names{1}; 
for i=2:length(f)
    pti_norm = pti./pt1; 
    if size(names)~=0
        pti = pt{i}; name = names{i}; 
        plot(st1(2:end).*TU, pti_norm(2:end), mkrs{i}, 'Color', 'blue', 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName', append(name, ' / ', name1));
    else
        plot(st1(2:end).*TU, pti_norm(2:end), mkrs{i}, 'Color', 'blue', 'LineWidth', 1, 'MarkerSize', 10, 'HandleVisibility', 'off');
    end
end
ylabel("Normalized program time", 'FontSize', 18, 'FontName', 'Times');

yyaxis right; hold on; 
s1 = s{1}; 
for i=2:length(f)
    si_norm = si./s1; 
    plot(st1(2:end).*TU, si_norm(2:end), mkrs{i}, 'Color', 'red', 'LineWidth', 1, 'MarkerSize', 10, 'HandleVisibility', 'off');
end
ylabel("Normalized cell count", 'FontSize', 18, 'FontName', 'Times');
xlim([0, st1(end).*TU]); 

ax.YAxis(1).Color = [0 0 1];
ax.YAxis(2).Color = [1 0 0];

title(tiles, SYS + " Runtime Comparison", 'FontWeight', 'Bold','FontName', 'Times','FontSize', 18);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [st, pt, s] = parse_time_files(f)
    st = {}; pt = {}; s = {}; 
    for i=1:length(f)
        sti = []; pti = sti; si = sti; fi = f{i};
        
        while ~feof(fi)
            line = split(fgetl(fi)); 
            pti(end+1) = str2double(line{5});
            sti(end+1) = str2double(line{9});
            size = strsplit(line{13},'/');
            si(end+1) = str2double(size{1}); 
        end
        
        st{end+1} = sti;
        pt{end+1} = pti;
        s{end+1}  = si;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%