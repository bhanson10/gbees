function plt = plot_nongaussian_surface(X,P,varargin)
% plot_nongaussian_surface.m
% Benjamin Hanson, 2024
% 
% Given a set of 2D/3D state vectors X with associated weights P,
% generate an isosurface representing a curve of isovalue
% 
% Inputs:
%          X -- set of 1D/2D/3D state coordinates
%          P -- weights of state vectors
%   varargin -- optional arguments
%               * p -- plotting parameters (optional)
%                   *   color -- isosurface color
%                   * display -- handle visibility
%                   *    name -- display name, if display==1
%                   *   means -- plot weighted mean of point mass PDF
%                   *     axh -- figure axis
%                   *   alpha -- surface visibility
%                   *    type -- distribution type
%               * aR -- alpha radius (for 2D and 3D functionality only)
%               * isovalue -- isosurface value(s) to plot (for 2D and 3D functionality only)
%               * axc -- figure axis to copy (because plotting alphaShape
%               automatically changes the aspect ratio and limits of a plot

% variable arguments - defaults
p.color = {"r", "r", "r"}; 
p.display = 0; 
p.means = 0;
p.axh = gca; 
p.alpha = [0.5, 0.3, 0.1];
p.type = "curve";
aR = [2 2 2]; 
isovalue = [0.68, 0.95, 0.997]; 
axc = p.axh; 

for i=1:2:length(varargin)
    if strcmp('p',varargin{i})
        p = varargin{i+1};
    elseif strcmp('aR',varargin{i})
        aR = varargin{i+1};
    elseif strcmp('isovalue',varargin{i})
        isovalue = varargin{i+1};
    elseif strcmp('axc',varargin{i})
        axc = varargin{i+1}; 
    else
        error(append("Unspecified argument: ", varargin{i}));
    end
end

% checks and balances
if length(X)~=length(P)
    error("Incongruous state vector/weight sets.")
end
if ~isfield(p,'color')
    for i = 1:numel(isovalue)
        p.color{i}="r";
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
    p.alpha=flip(logspace(log(0.5),log(0.75),numel(isovalue)));
else
    p.alpha = p.alpha.*ones(1,numel(isovalue));
end
if isfield(p,'type')
    if strcmp(p.type, "curve")
    elseif strcmp(p.type, "scatter")
        if isfield(p,'size')
            if isscalar(p.size)
                p.size = p.size.*ones(1,numel(isovalue));
            end
        else
            p.size = 10.*ones(1,numel(isovalue));
        end
    elseif strcmp(p.type, "grid")
        plt = plot_nongaussian_surface_g(X,P,'p',p,'axc',axc); 
        return; 
    else
        error("Unsupported type.")
    end
else
    p.type = "curve"; 
end
if exist("aR", "var")
    if isscalar(aR)
        if aR < 0
            error("alphaRadius may not be negative."); 
        else
            aR = aR.*ones(1,numel(isovalue));
        end
    else
        if numel(aR) ~= numel(isovalue)
            error("Incongruous alphaRadius/isovalue sets.");
        else
            if any(aR < 0)
                error("alphaRadius may not be negative."); 
            end
        end
    end
end
if (max(isovalue)>1)||(min(isovalue)<0)
    error("Isovalue is outside of probability bounds [0,1].")
end
if ~isa(axc, 'matlab.graphics.axis.Axes')
    error("Copy axis must be an axis variable.")
end

% Getting number of state vectors
N=size(X); 

switch N(2)
    case 1, plt = plot_nongaussian_surface1D(X,P,p);
    case 2, plt = plot_nongaussian_surface2D(X,P,p,aR,isovalue,axc);
    case 3, plt = plot_nongaussian_surface3D(X,P,p,aR,isovalue,axc);
    otherwise, error('Unsupported dimensionality');
end

function plt = plot_nongaussian_surface1D(X,P,p)
x_list = unique(X);
P_full=zeros(numel(x_list), 1);

for l=1:numel(P)
    i=find(x_list==X(l)); P_full(i)=P_full(i)+P(l);
end
P_full = P_full./sum(P_full); 

if p.display
    if strcmp(p.type, "curve")
        plt = plot(p.axh, x_list, P_full, "Color", p.color{1}, 'LineWidth', 2, 'DisplayName', p.name);
    elseif strcmp(p.type, "scatter")
        if strcmp(p.color{1}, "jet")
            plt = scatter(p.axh, x_list, P_full, p.size, P_full, 'filled', 'DisplayName', p.name);
            colormap(jet); 
            clim([min(P_full,[],'all'), max(P_full,[],'all')]);
        else
            plt = scatter(p.axh, x_list, P_full, p.size, p.color{1}, 'filled', 'DisplayName', p.name);
        end
    else
        error("Unsupported plot type.");
    end
else
    if strcmp(p.type, "curve")
        plt = plot(p.axh, x_list, P_full, "Color", p.color{1}, 'LineWidth', 2, 'HandleVisibility', 'off');
    elseif strcmp(p.type, "scatter")
        if strcmp(p.color{1}, "jet")
            plt = scatter(p.axh, x_list, P_full, p.size{1}, P_full, 'filled', 'HandleVisibility', 'off');
            colormap(jet); 
            clim([min(P_full,[],'all'), max(P_full,[],'all')]);
        else
            plt = scatter(p.axh, x_list, P_full, p.size{1}, p.color{1}, 'filled', 'HandleVisibility', 'off');
        end
    else
        error("Unsupported plot type.");
    end
end
if p.means
    x_mean = sum(P_full .* x_list); 
    if p.display
        xline(x_mean, "Color", p.color{1}, 'LineWidth', 2, 'DisplayName', p.name + " - mean");
    else
        xline(x_mean, "Color", p.color{1}, 'LineWidth', 2, 'HandleVisibility', "off");
    end
end

  
function plt = plot_nongaussian_surface2D(X,P,p,aR,isovalue,axc)
x_lim = axc.XLim;
y_lim = axc.YLim;
z_lim = axc.ZLim;
data_aspect_ratio = axc.DataAspectRatio; 
plot_aspect_ratio = axc.PlotBoxAspectRatio;

P_full = P ./ sum(P); 
[P_full, idx] = sort(P_full, 'descend');
x_list = X(idx, :);

if strcmp(p.type, "curve")
    Xs = {}; 
    for i = 1:numel(isovalue)
        Xs{i} = []; 
    end
    
    cum_sum = 0; idx = 1; 
    for i = 1:numel(P_full)
        if cum_sum > isovalue(idx)
            idx = idx + 1; 
            if idx > numel(isovalue)
                break
            end
        end
        for j = idx:numel(isovalue)
            Xs_j = Xs{j}; 
            Xs_j = [Xs_j; x_list(i, :)];
            Xs{j} = Xs_j; 
        end
        cum_sum = cum_sum + P_full(i); 
    end
end

for i = 1:numel(isovalue)
    if (i == 1)&&p.display
        if strcmp(p.type, "curve")
            Xi = Xs{i}; 
            if aR(i) > 0
                shp = alphaShape(Xi, aR(i));
            else
                shp = alphaShape(Xi);
            end
            plt{i} = plot(shp, "FaceColor", p.color{i}, "EdgeColor", "none", "FaceAlpha", p.alpha(i), "DisplayName", p.name, "Parent", p.axh);
        elseif strcmp(p.type, "scatter")
            if strcmp(p.color{1}, "jet")
                cmap = jet; bg = cmap(1, :); 
                set(p.axh, 'Color', bg);
                plt{i} = scatter(p.axh, x_list(:,1), x_list(:,2), p.size(i), P_full, 'filled', "DisplayName", p.name);
                colormap(jet); 
                clim([min(P_full,[],'all'), max(P_full,[],'all')]);
                colorbar;
            else
                plt{i} = scatter(p.axh, x_list(:,1), x_list(:,2), p.size(i), p.color{1}, 'filled', "DisplayName", p.name);
            end
            break;
        end
    else
        if strcmp(p.type, "curve")
            Xi = Xs{i}; 
            if aR(i) > 0
                shp = alphaShape(Xi, aR(i));
            else
                shp = alphaShape(Xi);
            end
            plt{i} = plot(shp, "FaceColor", p.color{i}, "EdgeColor", "none", "FaceAlpha", p.alpha(i), "HandleVisibility", "off", "Parent", p.axh);
        elseif strcmp(p.type, "scatter")
            if strcmp(p.color{1}, "jet")
                cmap = jet; bg = cmap(1, :); 
                set(p.axh, 'Color', bg);
                plt{i} = scatter(p.axh, x_list(:,1), x_list(:,2), p.size(i), P_full, 'filled', "HandleVisibility", "off");
                colormap(jet); 
                clim([min(P_full,[],'all'), max(P_full,[],'all')]);
                colorbar;
            else
                plt{i} = scatter(p.axh, x_list(:,1), x_list(:,2), p.size(i), p.color{1}, 'filled', "HandleVisibility", "off");
            end
            break;
        end
    end
end
if p.means, mean_X=sum(P .* x_list); scatter(p.axh, mean_X(:,1), mean_X(:,2), 100, p.color{1}, "pentagram", "filled", 'HandleVisibility', 'off'); end
if isequal(data_aspect_ratio, [1 1 1])
    axis equal
elseif isequal(plot_aspect_ratio, [1 1 1])
    axis square
else
    axis normal
end
p.axh.XLim = x_lim;
p.axh.YLim = y_lim;
p.axh.ZLim = z_lim;

function plt = plot_nongaussian_surface3D(X,P,p,aR,isovalue,axc)
x_lim = axc.XLim;
y_lim = axc.YLim;
z_lim = axc.ZLim;
data_aspect_ratio = axc.DataAspectRatio; 
plot_aspect_ratio = axc.PlotBoxAspectRatio;


P_full = P ./ sum(P); 
[P_full, idx] = sort(P_full, 'descend');
x_list = X(idx, :);

if strcmp(p.type, "curve")
    Xs = {}; 
    for i = 1:numel(isovalue)
        Xs{i} = []; 
    end
    
    cum_sum = 0; idx = 1; 
    for i = 1:numel(P_full)
        if cum_sum > isovalue(idx)
            idx = idx + 1; 
            if idx > numel(isovalue)
                break
            end
        end
        for j = idx:numel(isovalue)
            Xs_j = Xs{j}; 
            Xs_j = [Xs_j; x_list(i, :)];
            Xs{j} = Xs_j; 
        end
        cum_sum = cum_sum + P_full(i); 
    end
end

if p.means, mean_X=sum(P .* x_list); scatter3(p.axh, mean_X(:,1), mean_X(:,2), mean_X(:,3), 100, p.color{1}, "pentagram", "filled", 'HandleVisibility', 'off'); end

for i = 1:numel(isovalue)
    if (i == 1)&&p.display
        if strcmp(p.type, "curve")
            Xi = Xs{i}; 
            if aR(i) > 0
                shp = alphaShape(Xi, aR(i));
            else
                shp = alphaShape(Xi);
            end
            plt{i} = plot(shp, "FaceColor", p.color{i}, "EdgeColor", "none", "FaceAlpha", p.alpha(i), "DisplayName", p.name, 'Parent', p.axh);
        elseif strcmp(p.type, "scatter")
            if strcmp(p.color{1}, "jet")
                plt{i} = scatter3(p.axh, x_list(:,1), x_list(:,2), x_list(:,3), P_full, p.size(i), P_full, 'filled', "DisplayName", p.name);
                colormap(jet); 
                clim([min(P_full,[],'all'), max(P_full,[],'all')]);
            else
                plt{i} = scatter(p.axh, x_list(:,1), x_list(:,2), x_list(:,3), P_full, p.size(i), p.color{1}, 'filled', "DisplayName", p.name);
            end
        end
    else
        if strcmp(p.type, "curve")
            Xi = Xs{i}; 
            if aR(i) > 0
                shp = alphaShape(Xi, aR(i));
            else
                shp = alphaShape(Xi);
            end
            plt{i} = plot(shp, "FaceColor", p.color{i}, "EdgeColor", "none", "FaceAlpha", p.alpha(i), "HandleVisibility", "off", 'Parent', p.axh);
        elseif strcmp(p.type, "scatter")
            if strcmp(p.color{1}, "jet")
                plt{i} = scatter3(p.axh, x_list(:,1), x_list(:,2), x_list(:,3), P_full, p.size(i), P_full, 'filled', "HandleVisibility", "off");
                colormap(jet); 
                clim([min(P_full,[],'all'), max(P_full,[],'all')]);
            else
                plt{i} = scatter3(p.axh, x_list(:,1), x_list(:,2), x_list(:,3), P_full, p.size(i), p.color{1}, 'filled', "HandleVisibility", "off");
            end
        end
    end
end

if isequal(data_aspect_ratio, [1 1 1])
    axis equal
elseif isequal(plot_aspect_ratio, [1 1 1])
    axis square
else
    axis normal
end
p.axh.XLim = x_lim;
p.axh.YLim = y_lim;
p.axh.ZLim = z_lim;