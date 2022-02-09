% Run this script from within the directory containing the script file.

clear();

restoredefaultpath();
addpath(genpath("code"));

load("data\fig1_data.mat", "units", "tracking");
nUnits = numel(units);

%% Normalize and smooth firing rates, PCA

dt = 0.010;             % time bin size
speedSmoothSec = 1.0;   % gaussian smoothing for speed calculation (in seconds)
speedThresh = 0.025;    % exclude samples when animal is below this speed (in m/s)
sigmaFR = 5;            % gaussian smoothing of firing rate (in samples)

trngplt = 9114 + [0, 5]; % roger rec1 OF

% trng = session.timeRange;
validTimeRanges = [ 
    7457, 14778
    14890, 16045];

timeRangeAll = validTimeRanges([1, end]);
t = timeRangeAll(1) : dt : timeRangeAll(2);
t = t';
tEdges = cen2edg(t);
nt = numel(t);

[tplt, vplt] = restrict(t, trngplt);

% Calculate smoothed running speed of animal
trk = tracking;
sigmaSpeed = speedSmoothSec./dt;
nfilt = ceil(sigmaSpeed)*10 + 1;
xs = imgaussfilt(trk.x, sigmaSpeed, "filterSize", [nfilt, 1]);
ys = imgaussfilt(trk.y, sigmaSpeed, "filterSize", [nfilt, 1]);
trk.speed = [0; hypot(diff(xs), diff(ys)) ./ dt];
tracking = trk;

[~, trk.inValidTimeRange] = restrict(t, validTimeRanges);
trk.hasValidSpeed = trk.inValidTimeRange & trk.speed >= speedThresh;
trk.t = t;

v = trk.hasValidSpeed;
trk = structfun(@(x) x(v), trk, 'uni', 0);
t = t(v);
vplt = vplt(v);


X = zeros(nt, nUnits, 'single');
for u = 1:nUnits
    st = units(u).spikeTimes;
    X(:, u) = histcounts(st, tEdges);
end

X = imgaussfilt(X, sigmaFR, "filterSize", [51, 1]);

fprintf("Excluded %.1f%% bins below speed threshold\n", sum(~v)/nt*100);
X = X(v, :);

% Normalized firing rates
[Xz, Xmu, Xsigma] = zscore(X);
Xz = zscore(Xz, [], 2);
Xzplt = Xz(vplt, :);

% PCA
% Downsample the data points to every 25th sample (= 250 ms intervals)
subInterval = 5*sigmaFR;
isub = ceil( 1 : subInterval : sum(trk.hasValidSpeed) );
tsub = t(isub);
Xzsub = Xz(isub, :);
ntsub = numel(isub);

C = Xzsub'*Xzsub ./ (ntsub-1);

[V, d] = eig(C);
d = diag(d);
[d, isort] = sort(d, "descend");
V = V(:, isort);
explained = d ./ sum(d) * 100;
nkeep = 30;
explained = explained(1:nkeep);
V = V(:, 1:nkeep);

pcssub = Xzsub/V';
pcs = Xz/V';

clear X Xz


%% Plot all rate maps

clear allRateMaps

% figColScheme2 light
set(0, "defaultFigureWindowStyle", "normal");

fig = figure("position", [985 129 589 765], "colormap", viridis);

splx = 13;
sply = 15;

nbins = 50;
posRange = 1.50;

[~, yorder] = sort([units.spatialInfo], "descend");    % sort by spatial info
maxRates = zeros(nUnits, 1);

trk0 = tracking;
pos0 = double([trk0.x, trk0.y]);

[~, v] = restrict(trk0.t, validTimeRanges);
pos0 = pos0(v, :);
vspd = trk0.speed(v) > speedThresh;
pos0 = pos0(vspd, :);

t0 = trk0.t(v);
t0 = t0(vspd);

pc = prctile(pos0, [0.25, 99.75]);
lims = pc;

for n = 1:2
    edges{n} = linspace(lims(1, n), lims(2, n), nbins+1);
end

ndilate = 1;
rmsmooth = 2.75;
se = strel('disk', ndilate);

for u = 1:nUnits
    unit = units(yorder(u));
    st = unit.spikeTimes;
    st = restrict(st, validTimeRanges);
    sspd = interp1(trk0.t, trk0.speed, st);
    svalid = sspd > speedThresh;
    st = st(svalid);
    si = binarySearch(st, t0);
    rmf = ratemapq(pos0, si, edges, 0.01, "checkInputs", false, "smooth", rmsmooth, "ndilate", ndilate);
    alpha = rmf.validBin;
    
    ax = subplot(sply, splx, u, "parent", fig);
    resizeobj(ax, 1.3);
    imagesc(ax, rmf.z, "alphaData", alpha);
    ax.CLim = [0, prctile(rmf.z(:), 97.5)];
    maxRates(u) = max(rmf.z(:));
    allRateMaps(yorder(u)) = rmf;
end

cb = colorbar(ax, "location", "south");
cb.Units = "normalized";
cb.Position = [0.44353 0.095 0.21032 0.01];
cb.FontSize = 8;
ylabel(cb, "Normalized rate / A.U.");
cb.Ticks = ax.CLim;
cb.TickLabels = ["0", "max"];

axs = findobj(gcf, 'type', 'axes');
axis(axs, "image", "off")

%% Plot path in box

figure("position", [100 350 400 200]);
ax = gca;

x = 0.75 * [-1, -1, 1, 1, -1, -1];
y = 0.75 * [-1, 1, 1, -1, -1, 1];
plot(x, y, 'k', 'lineWidth', 2);
alpha = 0.05;

hold on
v = trk.inValidTimeRange;
plot(trk.x(v), trk.y(v), "color", [0, 0, 0, alpha], "lineWidth", 0.3);

v = vplt;
scatter(trk.x(v), trk.y(v), 5, 1:sum(v), "filled");
colormap(ax, "hot");
cb = colorbar(ax, "location", "eastoutside");
ylabel(cb, "Trajectory time / s");
cb.Ticks = cb.Limits;
cb.TickLabels = ["0", "5"];
cb.FontSize = 12;
cb.Position = [0.2 0.33 0.023 0.37];
axis(ax, "equal", "off");


%% Spike rasters and smoothed firing rates for all cells
% (Extended Data Fig. 3A

% Plot spike rates
figure("position", [788 183 400 600]);
sply = 7;
splx = 3;
icol = [1, 2];

vec = @(x) x(:);

isp = @(irow) vec((irow(:)-1)*splx + icol);

% Rasters
ax = subplot(sply, splx, isp([2, 3]));
x = [];
y = [];
for u = 1:nUnits
    st = units(yorder(u)).spikeTimes;
    st = restrict(st, trngplt);
    xtmp = st + [0, 0, nan];
    ytmp = zeros(size(st)) + [0, 1, nan] + u;
    x = [x; vec(xtmp')];
    y = [y; vec(ytmp')];
end
line(x, y, 'color', 'k', 'lineWidth', 0.5);
xticks([]);
ax.TickDir = "out";
ax.YTick = [1, nUnits];
axis tight
ax.YLim = [1, nUnits];
ylabel("Unit #");
hold on
ax.Clipping = 'off';

% colour 'bar' at top showing trajectory mapping
x = linspace(trngplt(1), trngplt(2), 1000);
y = zeros(size(x)) + nUnits + 10;
c = hot(numel(x));
scatter(x, y, 10, c, "filled");

% Smoothed rates
ax = subplot(sply, splx, isp([4, 5]));

vpltsub = vplt(isub);

x = tplt-tplt(1);
z = gather(Xzplt(:, yorder));
imagesc(tplt, 1:nUnits, z');
clim = prctile(z(:), [0.1, 99.9]);
ylabel("Unit #");
set(gca, 'clim', clim);
xticks([]);
ax.YDir = "normal";
ax.YTick = [1, nUnits];
ax.TickDir = "out";
colormap(ax, viridis());
cb = colorbar(ax, "location", "eastoutside", ...
    "position", [0.6898 0.42333 0.0252 0.0969]);
cb.Ticks = cb.Limits;
cb.TickLabels = ["min", "max"];
ylabel(cb, "Normalized rate / A.U.");

ax = subplot(sply, splx, isp(1));
hold on
ax.Color = [0, 0 ,0]+0.8;
y1 = gather(double(trk.x(vplt)));
y2 = gather(double(trk.y(vplt)));
plot(x, y1, 'color', 'k', 'lineWidth', 1);
plot(x, y2, 'color', [0.7, 0, 0], 'lineWidth', 1);
ylabel(ax, "Position / m");
ax.XTick = [];
text(ax, x(1)+0.2, y1(1)+0.05, "x", "color", "k", "horizontalAlignment", "left", "verticalAlignment", "bottom");
text(ax, x(1)+0.2, y2(1)+0.05, "y", "color", [0.7, 0, 0], "horizontalAlignment", "left", "verticalAlignment", "bottom");
ax.YLim = [-0.75, 0.75];
ax.YTick = [-0.75, 0, 0.75];

% PCs time course
ax = subplot(sply, splx, isp(6));
z = pcs(vplt, 1:6);
h = plot(x, z);
axis tight off
hleg = legend(h, string(1:6), "position", [0.69787     0.095278       0.1425      0.16917]);
title(hleg, "PC #");

%% Plot PC explained variance and spatial mapping

% Variance explained
if exist("explained", "var")
    fig = figure("colormap", viridis(),  "position", [652 200 300 300]);
    subplot(2, 2, 1);
    bar(explained, "faceColor", "k", "faceAlpha", 0.5, "edgeColor", "k");
    ylabel("% variance explained");
    xlabel("PC #");
    
    % Cumulative variance explained
    subplot(2, 2, 2);
    plot(cumsum(explained));
    ylabel("Cumulative % variance");
    xlabel("PC #");
    xline(6, "lineStyle", ":");
    
end

fig = figure("colormap", viridis(), "position", [741 323 460 440]);
pos = gather(double([t, trk.x, trk.y]));
sply = 4;
splx = 4;
nplt = sply*splx;

limits = [-1, 1, -1, 1]*0.75;

for n = 1:nplt
    subplot(sply, splx, n);
    resizeobj(gca, 1.2);
    z = [t, double(gather(pcs(:, n)))];
    z = gather(z);
    rm = analyses.map(pos, z, "limits", limits, "binWidth", 0.03, "smooth", 0);
    z = rm.z;
    z = smoothRateMap(z, rmsmooth);
    vbin = rm.countRaw > 0;
    vbin = imclose(vbin, se);
    imagesc(rm.x, rm.y, z, "alphaData", vbin);
    set(gca, 'clim', [-2.5, 2.5]);
end
axs = findall(fig, "type", "axes");
axis(axs, "off", "xy", "image");
cb = colorbar(axs(end), "position", [0.93 0.77 0.025 0.11]);
cb.YTick = sort([cb.Limits, 0]);

%% run UMAP on first 6 PCs

XU = double(gather(pcssub));
XU = XU(:, 1:6);

disp("Running UMAP. This will take a while ...");

umapArgs = {XU,  ...
    'n_components', 3, ...
    'n_neighbors', 5000, ...
    'metric', 'cosine', ...
    'randomize', false, ...
    'min_dist', 0.8, ...
    'repulsion_strength', 1, ...
    'method', 'c++', ...
    'parameter_names', cellstr(string(1:6)) };

[Xumap, umapObj] = run_umap(umapArgs{:});

disp("Finished running UMAP.");
disp("Attention!! Please note that the UMAP output rotation is arbitrary. This means you will need to rotate the 3-D plots to see the torus clearly!");


%% Scatter plot

showTraj = 1;
zoomfactor = 6;
camfov = 70;
markerSize = 10;

if showTraj
    % finely interpolate the subsampled UMAP trajectory for visualization
    Ptraj = Xumap(vpltsub, :);
    npsub = sum(vpltsub);
    ni = 5000;
    ii = linspace(1, npsub, ni)';
    Ptraji = interp1(Ptraj, ii, "pchip", nan);
    coltraj = hot(ceil(ni*1.2));
    coltraj(1:floor(ni*0.2), :) = [];
end

fig = figure();
fig.WindowStyle = "normal";
fig.Position = [200, 50, 512, 512];

P = gather(Xumap);
np = size(P, 1);

col = gather(pcssub(:, 1));
prc = prctile(col, [1, 99]);
col = interp1(linspace(prc(1), prc(2), 128), viridis(128), col, "linear", "extrap");

clear h
ax = gca;
resizeobj(ax, 1.5);
hold(ax, "on");
ax.Clipping = "off";
scatter3(ax, P(:, 1), P(:, 2), P(:, 3), markerSize, col, "filled", ...
    "markerFaceAlpha", 0.5, "markerEdgeColor", "none");

if showTraj
    htraj = scatter3(ax, Ptraji(:, 1), Ptraji(:, 2), Ptraji(:, 3), markerSize, coltraj, "filled", ...
        "markerFaceAlpha", 0.5, "markerEdgeColor", "none");
end

axis(ax, "equal", "off", "vis3d");
axis(ax, axis(ax));


%% Single-unit spikes on the torus

nex = 3;
maxspikes = 2000;
spikePosDither = 0.03;
scatterAlpha = 0.3;
scatterSzAll = 1.5;
scatterSzEx = 3;

fig = figure("position", [550, 300, 900, 280], "colormap", viridis);
clear axs
ex_id = ["2_1149", "2_1001", "2_1135"];

[v, inds] = ismember(ex_id, [units.id]);
assert(all(v));

for u = 1:nex
    unit = units(inds(u));
    st = unit.spikeTimes;
    st = restrict(st, timeRangeAll);
    sspd = interp1(trk0.t, trk0.speed, st, "linear", nan);
    svspd = sspd > 0.025;
    st = st(svspd);
    sp_torus = interp1(tsub, P, st, "linear", nan);
    sp_2d = interp1(trk0.t, [trk0.x, trk0.y], st, "linear", nan);
    if size(sp_torus, 1) > maxspikes
        irnd = randperm(size(sp_torus, 1), maxspikes);
        st = st(irnd);
        sp_torus = sp_torus(irnd, :);
        sp_2d = sp_2d(irnd, :);
    end
    nsp = numel(st);
    
    if spikePosDither
        sp_2d = sp_2d + 2*spikePosDither*(rand(nsp, 2)-0.5);
    end
    
    if u==1
        sp2d_of = sp_2d;
        spref_torus = sp_torus;
    end
    
    ax = subplot(1, 4, u);
    resizeobj(ax, 1.4);
    
    % plot all torus points
    h(u) = scatter3(P(:, 1), P(:, 2), P(:, 3), scatterSzAll, col, "filled", ...
        "markerEdgeColor", "none", "markerFaceAlpha", scatterAlpha);
    hold(ax, "on");
    
    % plot unit spike positions
    scatter3(ax, sp_torus(:, 1), sp_torus(:, 2), sp_torus(:, 3), scatterSzEx, "k", "filled", ...
        "markerEdgeColor", "none");
    
    text(0.95, 0.95, unit.id, 'sc', ...
        'horizontalAlignment', 'right', ...
        'verticalAlignment', 'top', ...
        'interpreter', 'none', ...
        'fontSize', 6);
    
    % INSET plots
    for n = 1:3
        
        % create inset axes
        axpos = ax.Position;
        leftpos = axpos(1) + (n-1)*0.06;
        axpos2 = [leftpos, axpos(2) + 0.1, 0.1, 0.18];
        ax2 = axes(fig, "position", axpos2);
        hold(ax2, "on");
        
        if n==1
            % spike+path plot
            ax2.XLim = limits(1:2);
            ax2.YLim = limits(3:4);
            plot(ax2, trk.x, trk.y, 'color', 0.8+[0, 0, 0]);
            plot(ax2, sp_2d(:, 1), sp_2d(:, 2), 'k.', 'markerSize', 2.5);
        elseif n==2
            % rate map
            rmap = allRateMaps(inds(u));
            z = rmap.z;
            v = rmap.validBin;
            z(~v) = nan;
            imagesc(ax2, z', "alphaData", v');
            ax2.CLim = [0, prctile(z(:), 99)];
            str = sprintf("%01.1f Hz", ax2.CLim(2));
            text(ax2, 0.5, 1.02, str, 'sc', "horizontalAlignment", "center", "verticalAlignment", "bottom", "fontSize", 6);
        elseif n==3
            % acorr
            rmap = allRateMaps(inds(u));
            v = rmap.validBin;
            z = rmap.z;
            z(~v) = mean(z(v));
            z = analyses.autocorrelation(z);
            gscore = analyses.gridnessScore(z);
            imagesc(ax2, z');
            ax2.CLim = [-1, 1];
            str = sprintf("GS = %01.2f", gscore);
            text(ax2, 0.5, 1.02, str, 'sc', "horizontalAlignment", "center", "verticalAlignment", "bottom", "fontSize", 6);
        end
        
        axis(ax2, "xy", "equal", "off");
        
    end
    
    
    axis(ax, "equal", "off", "vis3d");
    axs(u) = ax;
end

axis(axs, "equal", "off");