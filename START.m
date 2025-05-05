%%  Script performing density topology optimization for the minimization of characteristic number Lambda
% 
%  The method performs a local gradient-driven step evaluated by Adjoint
%  sensitivity analysis until it converges to a local minimum. The solution
%  space is regularized by the filtering technique (density and projection
%  filters) to penalize unwanted behavior. See [1] for a more detailed
%  description.
% 
%  The method uses matrices stemming from method-of-moments solution to field 
%  integral equation with Rao-Wilton-Glisson basis functions. 
%  In this script, the precalculated matrices result from 
%  electric field integral equation. Triangular mesh grid elements
%  are used as optimization degrees of freedom.
%
%  The state equation is Characteristic Modes analysis, a generalized
%  eigenvalue problem, which is evaluated in every iteration. The Adjoint
%  sensitivity analysis is developed based on properties of CMA, see [2].
% 
%  The magnitude of the first significant characteristic number is minimized on a design domain discretized into
%  800 triangles at the electrical size ka = 0.7
%
% [1] Tucek, J.,Capek, M., Jelinek, L., Sigmund, O.: Q-factor Minimization via Topology Optimization, 
%     pp. 1-13, 2023.
%
% [2] Tucek, J.,Capek, M., Jelinek, L.: Density-Based Topology Optimization for Characteristic Modes Manipulation, 
%     pp. 1-13, 2023.

clc; 
clear; 
close all;

% Precalculated matrices are loaded (they can be changed in AToM [1]):
load('operators_plate_ka0p7.mat');

% Auxiliarly matrices
OP.Z0    = Z0;      % vacuum part of the impedance matrix 
OP.Mesh  = Mesh;    % MATLAB structure of geometry (from AToM)
OP.BF    = BF;      % MATLAB structure of basis functions (from AToM)   
OP.BF2T  = BF2T;    % Connectivity matrix between triangles and basis functions

%% Fitness function settings

topOptSettings.nu = 0.1; % Parameter weighting the imaginary part of characteristic number
% Try nu = 1, nu = 0 and compare results.

% fitness function definition f = Re{Lambda}^2 + nu * Im{Lambda}^2
topOptSettings.fitness = @ff_minLam;

%% Topology optimization settings:
topOptSettings.interFun = @interFun;                  % Define interpolation function
topOptSettings.resistivityLimits = [1e5 0.01];        % Boundary resistivity values for vacuum (rho=0) and metal (rho=1)
topOptSettings.rmin = 0.15 * OP.Mesh.normDistanceA;   % Density filter radius
% The density filter is made larger to be more effective on a coarser mesh

% Projection filter is used in the continuation scheme
topOptSettings.beta = 1;       % Initial sharpness of the projection filter H
topOptSettings.etaVec  = 0.5;  % Level of the projection filter H
topOptSettings.betaIter = 75;  % Projection filter is updated every 75 iterations
topOptSettings.betaMax = 32;   % Maximal allowed sharpness of the projection filter H

topOptSettings.change = 0.01;  % Maximal change in rho between two consecutive runs
topOptSettings.maxIter = 6 * topOptSettings.betaIter;   % Maximal number of iterations

% Region near the feeder is fixed
topOptSettings.protTRs = [];    % Protected material (rho = 1)
topOptSettings.passiveTRs = []; % Passive material (rho = 0)

% Start the optimization
history = topOptWithCMA(OP, topOptSettings);

%% Postproccessing
% Convergence history
iter = 1:length(history.fval);
figure;
semilogy(iter, history.fval, 'x-', 'LineWidth', 1); hold on;
grid on;
xlabel('iteration');
ylabel('f (-)');
% add vertical dashed line when sharpness beta was updated
betaVec = unique(history.beta);
for i=1:length(betaVec)
    bb = find(betaVec(i) == history.beta);
    txt = ['it = ' num2str(bb(end)) ', beta = ' num2str(betaVec(i))];
    xline(bb(end),'--',{txt},'LabelHorizontalAlignment', 'left'); 
end
hold off;

legend('$f$',...
       '$\beta$ update','Interpreter','latex','FontSize',12)

% Plot final optimized design
iter = iter(end); % last iteration
betaPostProcess = history.beta(iter);
x = history.x(:,iter); % design variable 
xTilde = (history.H * x)./ sum(history.H,2);
xPhys = projectionFilter(xTilde, betaPostProcess, topOptSettings.etaVec); % physical (filtered) variable

plotDesign(OP.Mesh,xPhys);
[~,~,LamGray] = topOptSettings.fitness(OP, topOptSettings, xPhys);
title(['\nu = ' num2str(topOptSettings.nu,1) ', \lambda_1 = ' sprintf('%1.3f',real(LamGray)) '-j' sprintf('%1.3f',-imag(LamGray))])

% Hard thresholding of the gray design
xPhysBW = projectionFilter(xTilde, 150000, 0.5);

% Analyze BW design with physically accurate resistivites of vacuum and PEC
topOptSettings.resistivityLimits = [1e12 1e-8]; % Numerical values
[~,~,LamBW] = topOptSettings.fitness(OP, topOptSettings, xPhysBW);

plotDesign(OP.Mesh,xPhysBW); % plot structure
title(['\nu = ' num2str(topOptSettings.nu,1) ', \lambda_1 = ' sprintf('%1.3f',real(LamBW)) '-j' sprintf('%1.3f',-imag(LamBW))])
