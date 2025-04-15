function history = topOptWithCMA(OP, settings)
%% Perform the denstity topology optimization with Characteristic Modes Analysis
%   The magnitude of the first characteristic number is minimized
%   Refer to the paper [1] for more details about the density-based
%   topology optimization, and to [2] about its conjuction with
%   Characteristic Modes
% 
% Inputs:
%   OP           - MATLAB structure containing all variables and fields
%                  fully describing the optimization region and all
%                  necessary MoM matrices. Majority of the fields are
%                  evaluated by AToM, see [1]. Mandatory fields are:
%                  OP.Mesh (discretization), OP.BF (basis functions),
%                  OP.Z0 (vacuum impedance matrix), OP.BF2T (mapping matrix BF-TR)                 
%                  etc. See START.m for details.
%   settings     - MATLAB structure containing optimization settings.
%                  See START.m for details.
% 
% Outputs:
%   history      - MATLAB structure containing the history of the optimization process.
% 
% The code is initiated from START.m.
% 
% 2025, Jonas Tucek, CTU in Prague, jonas.tucek@fel.cvut.cz

%% Prepare Density filter
H = prepareDensityFilter(OP.Mesh,settings);

%% Initialize iteration
x = 0.5 * ones(OP.Mesh.nTriangles,1); % Initial design field distribution
x(settings.protTRs) = 1; % Protected triangles are metal
x(settings.passiveTRs) = 0; % Passive triangles are vacuum

% Indices of design variables (triangles)
settings.designTRs = setdiff(1:OP.Mesh.nTriangles,[settings.passiveTRs settings.protTRs]);

% Initial plot, which is to be updated
handle = plotDesign(OP.Mesh, x);

change = 1;
iter = 0;
loopbeta = 0;
nProtectedTRs = length(settings.protTRs);
nPassiveTRs = length(settings.passiveTRs);

%% Setup varibles for History struct
history.x = x;                  % Design variable
history.fval = [];              % fitness function values
history.change = [];            % change between consecutive iterations
history.beta = settings.beta;   % sharpness of the projection filter
history.Lam  = [];              % values of the first characteristic number
%% MMA setup
n = OP.Mesh.nTriangles ...
    - nProtectedTRs - nPassiveTRs;  % The number of design variables x_j.;
xmin  = zeros(n,1);                 % Column vector with the lower bounds for the variables x_j.
xmax  = ones(n,1);                  % Column vector with the upper bounds for the variables x_j.
xold1 = x;                          % xval, one iteration ago (provided that iter>1).
xold2 = x;                          % xval, two iterations ago (provided that iter>2).
low   = zeros(n,1);                 % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(n,1);                  % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).

AsymInc = 1.2;
AsymDecr = 0.7;

%% Start the topology optimization task

% Initial projection filter (with beta = 1 has zero influence)
[xPhys, dx] = projectionFilter(x, settings.beta, settings.etaVec);

tStart=tic;
while change > settings.change && iter < settings.maxIter
    iter = iter + 1;         % Iteration
    loopbeta = loopbeta + 1; % Number of iterations for projection filter update

    % Evaluate objectives and Adjoint sensitivity analysis
    [f, df, Lam] = settings.fitness(OP, settings, xPhys);

    % Use chain-rule to determine sensitivities with respect to the design variable x
    df = H * (df .* dx ./sum(H,2)); % Senstivities of objectives

    %% MMA update of the design variable
    % uncontrained version of MMA is utilized, see [2].
    xval = x(settings.designTRs);   % Only design triangles are updated
    df0dx = df(settings.designTRs); 

    % Perform MMA update
    [xmma, xold1, xold2, low, upp, AsymInc, AsymDecr, change] = ...
        mma_unconstrained(df0dx, xval, xold1, xold2, ...
        low, upp, iter, AsymInc, AsymDecr, xmin, xmax);
    
    % Retrieve updated design variable
    xnew = ones(OP.Mesh.nTriangles,1);
    xnew(settings.passiveTRs) = 0;
    xnew(settings.designTRs) = xmma;
    
    % Perform filtering
    xTilde = (H * xnew) ./sum(H,2); % Density filter
    [xPhys, dx] = projectionFilter(xTilde, settings.beta, settings.etaVec); % Projection filter

    x = xnew;  % New design variable

    %% Print and plot
    TrianglesColour = ([255 255 255] - xPhys*([255 255 255] - [0 0 0]))/255;
    handle.FaceVertexCData = TrianglesColour;
    pause(0.01);
    fprintf('iter.:%5i change:%7.3f f:%1.3f, Lam:%1.3f-j*%1.3f:\n',iter, ...
         change, f, real(Lam), -imag(Lam));

    %% Save History
    history.beta = [history.beta settings.beta];
    history.x = [history.x x];
    history.fval = [history.fval f];
    history.Lam = [history.Lam Lam];
    history.change = [history.change change];

    %% Update sharpness of the projection filter
    if settings.beta < settings.betaMax && (loopbeta >= settings.betaIter || change <= settings.change)
        settings.beta = 2 * settings.beta; % Double sharpness value
        change = 1;
        loopbeta = 0;
        fprintf('Sharpness of the projection filter, beta, increased to %g.\n',settings.beta);
    end
end
t0 = toc(tStart);
fprintf('\nTopology optimization found the design with Lambda=%1.3f in %3i iterations and %1.2f seconds\n',...
       Lam, iter, t0); 

history.t0 = t0;
history.H = H;
