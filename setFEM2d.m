function fem = setFEM2d(n)
% FEM = SETUPFEM2D(NEL)
% setup 2D FEM structure (uniform triangular grid)
% Input:  n   -- number of grid points in 1D
% Output: fem -- structure containing
%          nel   -- number of elements
%          xx,yy -- nodes
%          h2    -- mesh constant: 2*area(T)
%          A     -- function A(u) returning stiffness matrix (-\Delta + u)
%          M     -- mass matrix for right hand side in V_h
%          Y     -- function Y(y) returning weighted mass matrix
%          P     -- projection V_h to U_h: trapezoidal rule
%          obs   -- nodes of observation boundary
%
% Christian Clason (christian.clason@uni-graz.at)
% Bangti Jin       (btjin@math.tamu.edu)
% February 27, 2011

%% Set up grid: uniform triangular mesh
x       = linspace(0,1,n)';      % spatial grid points (uniform in x,y)
h       = (x(2)-x(1)); h2 = h^2; % Jacobi determinant of transformation (2*area(T))
xc      = x(2:end)-h/2;          % cell centers
[xx,yy] = meshgrid(x);           % coordinates of nodes
nel     = 2*(n-1)^2;             % number of nodes
tri     = zeros(nel,3);          % triangulation
ind = 1;
for i = 1:n-1
    for j = 1:n-1
        node         = (i-1)*n+j+1;              % two triangles meeting at node
        tri(ind,:)   = [node node-1 node+n];     % first triangle (lower left)
        tri(ind+1,:) = [node+n-1 node+n node-1]; % second triangle (upper right)
        ind = ind+2;
    end
end
% Creating rhs matrix

%ubdnod = find(xx == 0);     % unknown boundary nodes flipud(find(xx(:) == 0))
ubdnod = [flipud(find(xx(:) == 1))];
ybdnod = [flipud(find(yy(:) == 0)) flipud(find(yy(:) == 1)) flipud(find(xx(:) == 0))]; % observation boundary nodes

%% stiffness and mass matrix
 [K,L] = assembleFEM(tri,h2);
     M = assembleM(ubdnod,h,n^2);

%% projection on piecewise constants
%P0    = @(y) (y(ubdnod(1:end-1))+y(ubdnod(2:end)))/2;

%% FEM structure: grid, stiffness & mass matrices, projections etc.
fem.nel = nel;                         % number of elements
fem.xx  = xx;                          % grid
fem.yy  = yy;                          % grid
fem.tri = tri;
fem.xc  = xc;                          % cell centers on Robin boundary
fem.h2  = h2;                          % mesh constant: 2*area(T)
fem.A   = K + assembleS(ybdnod,h,n^2); % differential operator
fem.M   = M;                           % mass matrix on boundary
fem.L   = L;                           % mass matrix in domain
%fem.Y   = assembleY(ubdnod,h,n^2);     % weighted mass matrix: Y(y)*du = <y du,v>
%fem.P   = P0;                          % projection V_h to U_h
fem.obs = ybdnod(:);                   % trace on observation boundary
end % function setupFEM2D


%% Assemble stiffness matrix
function K = assembleK(t)
nel = size(t,1);    % number of elements

Ke = 1/2*[2 -1 -1 -1 1 0 -1 0 1]'; % elemental stiffness matrix <phi_i',phi_j'>

ent = 9*nel;
row = zeros(ent,1);
col = zeros(ent,1);
valk = zeros(ent,1);

ind=1;
for el=1:nel
    ll      = ind:(ind+8);         % local node indices
    gl      = t(el,:);             % global node indices
    cg      = gl([1;1;1],:); rg = gl';
    rg      = rg(:,[1 1 1]);
    row(ll) = rg(:);
    col(ll) = cg(:);
    valk(ll)= Ke;
    ind     = ind+9;
end
K = sparse(row,col,valk);
end % function assembleK

%% Assemble mass matrix <u,v>_Gamma_c
function U = assembleM(bdnod,h,n2)
nel = size(bdnod,1)-1;             % number of elements on boundary
nbd = size(bdnod,2);               % number of boundaries
Me = h/6 * [2 1 1 2]';             % elemental mass matrix <phi_i,phi_j>
%
ent = 4*nel*nbd;
row = zeros(ent,1);
col = zeros(ent,1);
val = zeros(ent,1);

ind=1;
for bd=1:nbd;                      % integrate over each subboundary
    for el=1:nel
        gl      = ind:ind+3;
        row(gl) = [bdnod(el,bd) bdnod(el+1,bd) bdnod(el,bd) bdnod(el+1,bd)];
        col(gl) = [bdnod(el,bd) bdnod(el,bd) bdnod(el+1,bd) bdnod(el+1,bd)];
        val(gl) = Me;
        ind     = ind+4;
    end
end
U = sparse(row,col,val,n2,n2);
end % function assembleM

%% Assemble potential mass matrix <y,v>_Gamma_i for u in U_h
function U = assembleS(bdnod,h,n2)
nel = size(bdnod,1)-1;             % number of elements on boundary

Me = h/6 * [2 1 1 2]';             % elemental mass matrix <phi_i,phi_j>

ent = 4*nel;
row = zeros(ent,1);
col = zeros(ent,1);
val = zeros(ent,1);

ind=1;
for el=1:nel
    gl      = ind:ind+3;
    row(gl) = [bdnod(el) bdnod(el+1) bdnod(el) bdnod(el+1)];
    col(gl) = [bdnod(el) bdnod(el) bdnod(el+1) bdnod(el+1)];
    val(gl) = Me;
    ind     = ind+4;
end
U = sparse(row,col,val,n2,n2);
end % function assembleU


%% Assemble weighted mass matrix: Y(y)*du = <y du,v> for du in U_h
%% Assemble stiffness, mass matrix
function [K,L] = assembleFEM(t,h2)
nel = size(t,1);    % number of elements

Ke = 1/2*[2 -1 -1 -1 1 0 -1 0 1]'; % elemental stiffness matrix <phi_i',phi_j'>
Me = h2/24 * [2 1 1 1 2 1 1 1 2]'; % elemental mass matrix <phi_i,phi_j>

ent = 9*nel;
row = zeros(ent,1);
col = zeros(ent,1);
valk = zeros(ent,1);
valm = zeros(ent,1);

ind=1;
for el=1:nel
    ll      = ind:(ind+8);         % local node indices
    gl      = t(el,:);             % global node indices
    cg      = gl([1;1;1],:); rg = gl';
    rg      = rg(:,[1 1 1]);
    row(ll) = rg(:);
    col(ll) = cg(:);
    valk(ll)= Ke;
    valm(ll)= Me;
    ind     = ind+9;
end
K = sparse(row,col,valk);
L = sparse(row,col,valm);
end % function assembleFEM
