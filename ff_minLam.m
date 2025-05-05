function [f, df, Lam] = ff_minLam(OP, settings, x)
%% Evalute the first characteristic number and compute sensitivities by Adjoint Sensitivity analysis
%       This function evaluates lambda_1 and punishes its imaginary part
%
% Inputs:
%   OP              - MATLAB structure containing all required information
%                     about the structure. Matrices Z0, and BF2T are required.
%   settings        - MATLAB structure containing optimization settings.
%                     Namely, boundary resistivities, normalization constant and design triangles labels, is necessary.
%   x               - design variable [0,1]^Tx1
%
% Outputs:
%   f      - fitness function value
%   df     - derivative of the fitness function
%   Lam    - value of the first characteristic number 
%
% 2025, Jonas Tucek, CTU in Prague, jonas.tucek@fel.cvut.cz

%% MoM Analysis
[Zs, dZs] = settings.interFun(x, settings.resistivityLimits); % Map design variable to resistivities
Zrho = ConstructMaterialMatrix(OP.Mesh, OP.BF, Zs, 1:OP.Mesh.nTriangles); % Construct the material matrix Zrho

% Evaluate Characteristic Modes
R = real(OP.Z0);
X = imag(OP.Z0) - 1j*Zrho;

[I, Lam] = eigs(X, R, 1, 'sm');
I = I./sqrt(I'*R*I); % Normalize

% Evaluate the objective
f = real(Lam)^2 + settings.nu * imag(Lam)^2;

if nargout > 1 % Adjoint sensitivity analysis provided on-fly

    LHS = [(X-Lam*R), -R*conj(I);...
        -I.'*R, 0];
    RHS = [zeros(OP.BF.nUnknowns,1);...
        -1/2 * (2 * real(Lam) - 1j * 2 * settings.nu * imag(Lam))];

    adj = LHS \ RHS; % compute adjoint variable

    z = adj(1:end-1); % to evaluate sensitivities only "z" adjoint is required

    df = zeros(OP.Mesh.nTriangles, 1);
    % Evaluate sensitivities
    for iT = 1:OP.Mesh.nTriangles  % Sweep through design triangles
        BFs = find(OP.BF2T(iT,:)); % find corresponding basis functions of the triangle
        Ie = I(BFs);               % eigencurrent flowing through corresponding basis functions
        z_e = z(BFs);

        dZrhoEl = ConstructMaterialMatrix(OP.Mesh,OP.BF, dZs, iT); % Compute element matrix
        % Final sensitivities
        df(iT,:) = 0 - 2*real(1j*z_e.' * dZrhoEl(BFs,BFs) * Ie);
    end
end


