function plot_nongaussian_surface(X,P,isovalue,p)
% plot_nongaussian_surface.m
% Benjamin Hanson, 2024
% 
% Given a set of 2D/3D state vectors X with associated weights P,
% generate an isosurface representing a curve of isovalue
% 
% Inputs:
%          X -- set of 2D/3D state vectors
%          P -- weights of 2D/3D state vectors
%   isovalue -- isosurface value(s) to plot
%          p -- plotting parameters (optional)
%               *   color -- isosurface color
%               * display -- handle visibility
%               *    name -- display name, if display==1
%               *   means -- plot weighted mean of point mass PDF
%               *     axh -- figure axis
%               *   alpha -- surface visibility
%               *    type -- distribution type
% plot_nongaussian_surface.m, https://github.com/bhanson10/gbees/tree/main/examples
% Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

% Checks and Balances
if length(X)~=length(P)
    error("Incongruous state vector/weight sets.")
end
if (max(isovalue)>1)||(min(isovalue)<0)
    error("Isovalue is outside of probability bounds [0,1].")
end
if ~exist('p','var')
    for i = 1:numel(isovalue)
        p.color{i}=[1 0 0];
    end
    p.display=0; 
    p.means=0;
    p.axh=gca; 
    p.alpha=flip(logspace(log(0.3),log(0.6),numel(isovalue)));
    p.type = "grid";
else
    if ~isfield(p,'color')
        for i = 1:numel(isovalue)
            p.color{i}=[1 0 0];
        end
    else
        if (isstring(p.color))||(ischar(p.color))||((all(size(p.color) == [1,3]))&&(~iscell(p.color)))
            col = p.color; p.color = {}; 
            for i = 1:numel(isovalue)
                p.color{i}=col;
            end
        end 
    end

    if ~isfield(p,'display')
        p.display=0;
    end
    if(p.display == 1)
        if ~isfield(p,'name')
            p.display=0;
        end
    end
    if isfield(p,'name')
        if ~isfield(p,'display')
            p.display=1;
        end
    end
    if ~isfield(p,'means')
        p.means=0;
    end
    if ~isfield(p,'axh')
        p.axh=gca;
    end
    if ~isfield(p,'alpha')
        p.alpha=flip(logspace(log(0.3),log(0.6),numel(isovalue)));
    else
        p.alpha = p.alpha.*ones(1,numel(isovalue));
    end

    if ~isfield(P, "type")
        p.type = "grid";
    end
end
p.alpha = sort(p.alpha); 
isovalue = sort(isovalue); 

% Getting number of state vectors
N=size(X); 

switch N(2)
    case 2, plot_nongaussian_surface2D(X,P,isovalue,p);
    case 3, plot_nongaussian_surface3D(X,P,isovalue,p);
   otherwise
      error('Unsupported dimensionality');
end

function plot_nongaussian_surface2D(X,P,isovalue,p)

if strcmp(p.type, "grid")
    x_list = unique(X(:,1)); x_list = [2*x_list(1) - x_list(2); x_list; 2*x_list(end) - x_list(end-1)];
    y_list = unique(X(:,2)); y_list = [2*y_list(1) - y_list(2); y_list; 2*y_list(end) - y_list(end-1)];
    [X_grid,Y_grid] = meshgrid(x_list, y_list);
    P_full=zeros(numel(y_list), numel(x_list));
    
    for l=1:numel(P)
        i=find(x_list==X(l,1)); j=find(y_list==X(l,2)); 
        P_full(j,i)=P_full(j,i)+P(l);
    end
    P_full = P_full./max(P_full,[],'all'); 
    
    count = 1; 
    for i=isovalue
        if (count == 1)&&p.display
            contour(p.axh, X_grid, Y_grid, P_full, [i i], '-', 'EdgeAlpha', p.alpha(count), 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count), 'FaceColor', p.color{count}, 'Fill', 'on', 'DisplayName', p.name);
        else
            contour(p.axh, X_grid, Y_grid, P_full, [i i], '-', 'EdgeAlpha', p.alpha(count), 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count), 'FaceColor', p.color{count}, 'Fill', 'on', 'HandleVisibility', 'off');
        end
        count = count + 1; 
    end
end

if p.means, mean_X=X.*P; scatter(p.axh, mean(mean_X(:,1)), mean(mean_X(:,2)), 100, p.color{1}, "pentagram", "filled", 'HandleVisibility', 'off'); end

function plot_nongaussian_surface3D(X,P,isovalue,p)

if p.means, mean_X=sum(X.*P); scatter3(p.axh, mean_X(1), mean_X(2), mean_X(3), 100, p.color{1}, "pentagram", "filled", 'HandleVisibility', 'off'); end

if strcmp(p.type, "grid")
    x_list = unique(X(:,1)); x_list = [2*x_list(1) - x_list(2); x_list; 2*x_list(end) - x_list(end-1)];
    y_list = unique(X(:,2)); y_list = [2*y_list(1) - y_list(2); y_list; 2*y_list(end) - y_list(end-1)];
    z_list = unique(X(:,3)); z_list = [2*z_list(1) - z_list(2); z_list; 2*z_list(end) - z_list(end-1)]; 
    [X_grid,Y_grid,Z_grid] = meshgrid(x_list, y_list, z_list);
    P_full=zeros(numel(y_list), numel(x_list), numel(z_list));
    
    for l=1:numel(P)
        i=find(x_list==X(l,1)); j=find(y_list==X(l,2)); k=find(z_list==X(l,3)); 
        P_full(j,i,k)=P_full(j,i,k)+P(l);
    end
    P_full = P_full./max(P_full,[],'all');
    
    count = 1; 
    for i=isovalue
        if (count == 1)&&p.display
            patch(isosurface(X_grid, Y_grid, Z_grid, P_full, i), 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count), 'FaceColor', p.color{count},'DisplayName', p.name); 
        else
            patch(isosurface(X_grid, Y_grid, Z_grid, P_full, i), 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count), 'FaceColor', p.color{count},'HandleVisibility', 'off');
        end
        count = count + 1; 
    end
end