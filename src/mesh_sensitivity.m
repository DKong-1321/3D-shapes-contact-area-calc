% meshSensitivity_generateMesh_autoStop.m
% Runs mesh refinement until the % change in contact area between consecutive
% meshes is < 5% for 3 consecutive refinements (i.e., once it first drops
% below 5%, it will still run 2 more finer meshes). Also saves a contact
% visualisation for each mesh (sphere mesh + highlighted contact triangles).



% PATHS
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
modelDir = fullfile(projectRoot,'model');

planeStlPath  = fullfile(modelDir,'50mm cube 4 iterations.stl');
sphereStlPath = fullfile(modelDir,'50mm diameter sphere.stl');

figDir = scriptDir;

% SETTINGS
opts = struct();
opts.tol = 0.20;
opts.sampleThreshold   = 0.5;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 30;
opts.roiExpandFactor   = 1.5;

errThreshPct = 5;
needConsec   = 3;

Hmax0        = 5.00;
refineFactor = 0.80;
HmaxMin      = 0.10;
maxIters     = 30;

makeFiguresPerMesh = true;

saveFigPerMesh = true;   % saves .fig
savePngPerMesh = true;   % saves .png

showFigPerMesh = true;   % SHOW LIVE FIGURES WHILE RUNNING
pausePerMesh   = 0.15;   % seconds (set 0 to avoid pauses)

figSphereAlpha     = 0.25;
figSphereColor     = [0 0.7 1];
figContactColor    = [0 1 0];
figMasterColor     = [0.85 0.85 0.85];
figEdgeColor       = [0.2 0.2 0.2];

placeSphereOnTop = true;
zGap = 0.0;

% LOAD MASTER (FIXED)
[F_plane, V_plane] = loadStlMesh(planeStlPath, 'plane');
master = buildBodyStruct(F_plane, V_plane);

% STORAGE (preallocate)
HmaxVals    = nan(maxIters,1);
sphereFaces = nan(maxIters,1);
meanEdge    = nan(maxIters,1);
Acontact    = nan(maxIters,1);
dA_pct      = nan(maxIters,1);

consecGood = 0;
iter = 0;
Hmax = Hmax0;

k = 0;
Aprev = NaN;

% LOOP UNTIL 3 CONSECUTIVE < 5%
while true
    iter = iter + 1;
    if iter > maxIters
        warning('Reached maxIters without satisfying convergence criterion.');
        break;
    end
    if Hmax < HmaxMin
        warning('Reached HmaxMin without satisfying convergence criterion.');
        break;
    end

    fprintf('Mesh level %d | Hmax = %.4f mm\n', iter, Hmax);

    model = femodel(Geometry=sphereStlPath);
    model = generateMesh(model, "Hmax", Hmax, "GeometricOrder", "linear");

    msh = model.Geometry.Mesh;
    TR  = triangulation(msh.Elements', msh.Nodes');
    [F_surf, V_surf] = freeBoundary(TR);
    [F_surf, V_surf] = patchCleanUnused(F_surf, V_surf);

    if placeSphereOnTop
        V_surf = placeMeshOnTopOfMaster(master.V, V_surf, zGap);
    end

    stats = meshEdgeLengthStats(F_surf, V_surf);

    slave = buildBodyStruct(F_surf, V_surf);
    results = computeContactArea_STS(master, slave, opts);

    k = k + 1;
    HmaxVals(k)    = Hmax;
    sphereFaces(k) = size(F_surf,1);
    meanEdge(k)    = stats.mean;
    Acontact(k)    = results.contactArea;

    if k == 1
        dA_pct(k) = NaN;
    else
        dA_pct(k) = 100 * abs(Acontact(k) - Aprev) / max(eps, abs(Aprev));
    end
    Aprev = Acontact(k);

    if k >= 2 && dA_pct(k) < errThreshPct
        consecGood = consecGood + 1;
    else
        consecGood = 0;
    end

    if makeFiguresPerMesh
        contactIdx = extractSlaveContactTriIdx(results, size(F_surf,1));
        saveContactFigure(figDir, k, Hmax, master, F_surf, V_surf, contactIdx, results.contactArea, ...
            figMasterColor, figSphereColor, figContactColor, figEdgeColor, figSphereAlpha, ...
            saveFigPerMesh, savePngPerMesh, showFigPerMesh, pausePerMesh);
    end

    if consecGood >= needConsec
        fprintf('Converged: ΔA < %.1f%% for %d consecutive refinements.\n', errThreshPct, needConsec);
        break;
    end

    Hmax = Hmax * refineFactor;
end


HmaxVals    = HmaxVals(1:k);
sphereFaces = sphereFaces(1:k);
meanEdge    = meanEdge(1:k);
Acontact    = Acontact(1:k);
dA_pct      = dA_pct(1:k);

% RESULTS TABLE
T = table(HmaxVals, sphereFaces, meanEdge, Acontact, dA_pct, ...
    'VariableNames', {'Hmax_mm','SphereFaces','MeanEdge_mm','ContactArea_mm2','DeltaArea_pct'});
disp(T);

% PLOTS
figure;
semilogx(sphereFaces, Acontact, '-o', 'LineWidth', 1.5);
xlabel('Number of sphere faces (log scale)');
ylabel('Contact area (mm^2)');
grid on;
exportgraphics(gcf, fullfile(figDir, 'meshSensitivity_contactArea_vs_numFaces.png'), 'Resolution', 300);

figure;
semilogx(sphereFaces(2:end), dA_pct(2:end), '-o', 'LineWidth', 1.5);
yline(errThreshPct, '--');
xlabel('Number of sphere faces (log scale)');
ylabel('\Delta Contact area (%)');
grid on;
exportgraphics(gcf, fullfile(figDir, 'meshSensitivity_deltaArea_vs_numFaces.png'), 'Resolution', 300);

% HELPERS
function [F2,V2] = patchCleanUnused(F,V)
    used = unique(F(:));
    map  = zeros(max(used),1);
    map(used) = 1:numel(used);
    V2 = V(used,:);
    F2 = map(F);
end

function stats = meshEdgeLengthStats(F,V)
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E,2);
    E = unique(E,'rows');
    d = V(E(:,1),:) - V(E(:,2),:);
    L = sqrt(sum(d.^2,2));
    stats.mean = mean(L);
end

function V_slave2 = placeMeshOnTopOfMaster(V_master, V_slave, zGap)
    mnM = min(V_master,[],1); mxM = max(V_master,[],1);
    mnS = min(V_slave ,[],1); mxS = max(V_slave ,[],1);

    cM = 0.5*(mnM + mxM);
    cS = 0.5*(mnS + mxS);

    dx = cM(1) - cS(1);
    dy = cM(2) - cS(2);
    dz = mxM(3) - mnS(3) + zGap;

    V_slave2 = V_slave + [dx dy dz];
end

function idx = extractSlaveContactTriIdx(results, nFacesSlave)
    idx = [];
    if isstruct(results)
        cand = {};
        if isfield(results,'contactTriIdxSlave'), cand{end+1} = results.contactTriIdxSlave; end
        if isfield(results,'slaveContactTriIdx'), cand{end+1} = results.slaveContactTriIdx; end
        if isfield(results,'slaveContactIdx'),    cand{end+1} = results.slaveContactIdx;    end
        if isfield(results,'contactTriIdx'),      cand{end+1} = results.contactTriIdx;      end
        if isfield(results,'contactFacesSlave'),  cand{end+1} = results.contactFacesSlave;  end
        if isfield(results,'slaveContactFaces'),  cand{end+1} = results.slaveContactFaces;  end

        for i = 1:numel(cand)
            x = cand{i};
            if isempty(x), continue; end

            if islogical(x)
                if numel(x) == nFacesSlave, idx = find(x); return; end
            end

            if isnumeric(x)
                x = x(:);
                x = x(isfinite(x));
                x = unique(round(x));
                x = x(x >= 1 & x <= nFacesSlave);
                if ~isempty(x), idx = x; return; end
            end
        end
    end
end

function saveContactFigure(figDir, k, Hmax, master, F_s, V_s, contactIdx, A, ...
    masterColor, sphereColor, contactColor, edgeColor, sphereAlpha, saveFig, savePng, showLive, pauseSeconds)

    if showLive
        f = figure('Color','w');
    else
        f = figure('Visible','off','Color','w');
    end

    ax = axes(f); hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
    view(ax, 3);

    patch('Faces', master.F, 'Vertices', master.V, ...
        'FaceColor', masterColor, 'EdgeColor', 'none', 'FaceAlpha', 1.0);

    patch('Faces', F_s, 'Vertices', V_s, ...
        'FaceColor', sphereColor, 'EdgeColor', edgeColor, ...
        'FaceAlpha', sphereAlpha, 'LineWidth', 0.25);

    if ~isempty(contactIdx)
        patch('Faces', F_s(contactIdx,:), 'Vertices', V_s, ...
            'FaceColor', contactColor, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
    end

    title(ax, sprintf('Mesh %d | Hmax=%.3f mm | A=%.3f mm^2 | contact faces=%d', ...
        k, Hmax, A, numel(contactIdx)), 'Interpreter','none');
    xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
    camlight(ax,'headlight'); lighting(ax,'gouraud');

    rotate3d(f,'on');

    baseName = sprintf('mesh_%02d_Hmax_%0.3f_contact', k, Hmax);

    if saveFig
        try
            savefig(f, fullfile(figDir, [baseName '.fig']));
        catch ME
            warning('Could not save .fig: %s', ME.message);
        end
    end
    if savePng
        try
            exportgraphics(f, fullfile(figDir, [baseName '.png']), 'Resolution', 250);
        catch ME
            warning('Could not save .png: %s', ME.message);
        end
    end

    if showLive
        drawnow;
        if pauseSeconds > 0
            pause(pauseSeconds);
        end
    else
        close(f);
    end
end
