function handle = plotDesign(Mesh, x)
%% Function plotting the optimized design
%
% Inputs:
%   Mesh        - MATLAB structure of geometry (from AToM, see [1])
%   x           - design variable [0,1]^Tx1
%
% Outputs:
%   handle      - handle to the graphical object, which is used to update
%                 the plot
%
% 2025, Jonas Tucek, CTU in Prague, jonas.tucek@fel.cvut.cz

% design variable is mapped to color [white, black] = [0, 1]
vecTrianglesColour = ([255 255 255] - x*([255 255 255] - [0 0 0]))/255;

figure('color', 'w');
handle=trisurf(Mesh.connectivityList, ...
    Mesh.nodes(:,1), Mesh.nodes(:,2), Mesh.nodes(:,3));
view(0, 90);
axis equal;
xlim([min(Mesh.nodes(:,1)) max(Mesh.nodes(:,1))]);
ylim([min(Mesh.nodes(:,2)) max(Mesh.nodes(:,2))]);

handle.FaceColor = 'flat';
handle.CDataMapping = 'scaled';
handle.FaceVertexCData = vecTrianglesColour;

ax = gca;
ax.XTick = [];
ax.YTick = [];
end
