clear all; close all; clc; 

mkrs = {"-o", "-square", "->", "-pentagram"}; 
MDL = "<model_type>";
DIR = {"c", "python"}; 
TU = 1; 
names = {'c implementation','python wrapper'};
time_unit = '(TU)';

f = {}; 
for i=1:length(names)
    f{end+1} = fopen("./" + MDL + "/results/" + DIR{i} + "/runtime.txt", 'r');
end

[st, pt, s] = parse_time_files(f);

f1 = figure(1); clf; f1.Position = [100 350 1300 475];
tiles = tiledlayout(1,2);

nexttile; ax = gca; 
l = legend; l.Location = "Northwest"; l.FontSize = 14; l.FontName = "Times"; 
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
ylabel("Number of cells", 'FontSize', 18, 'FontName', 'Times');
xlim([0, sti(end).*TU]); 

nexttile; ax = gca; 
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
    if size(names)~=0
        si = s{i}; name = names{i};  
        plot(st1(2:end).*TU, si_norm(2:end), mkrs{i}, 'Color', 'red', 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName', append(name, ' / ', name1));
    else
        plot(st1(2:end).*TU, si_norm(2:end), mkrs{i}, 'Color', 'red', 'LineWidth', 1, 'MarkerSize', 10, 'HandleVisibility', 'off');
    end
end
ylabel("Normalized cell number", 'FontSize', 18, 'FontName', 'Times');
xlim([0, st1(end).*TU]); 

ax.YAxis(1).Color = [0 0 1];
ax.YAxis(2).Color = [1 0 0];

title(tiles, MDL + " Runtime Comparison", 'FontWeight', 'Bold','FontName', 'Times','FontSize', 18);
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