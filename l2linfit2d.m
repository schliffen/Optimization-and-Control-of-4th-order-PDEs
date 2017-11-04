function u = l2linfit2d(fem,Mf,yd,alpha)
% [U,P] = L2NONLINFIT(FEM,MF,YD,ALPHA,UE,OPTS)
% solve nonlinear L2 data fitting problem using Newton method
%
% Input:  fem    -- finite element structure containing matrices
%                   (A(u),M,Y) and projection operators (P,P2)
%            Mf     -- right hand side with mass matrix already applied
%            yd     -- measured data
%            alpha  -- regularization parameter for L2 penalty
%            ue     -- exact solution (to evaluate reconstruction error)
%            opts   -- options structure containing (defaults)
%            maxit -- max. number of Newton iterations (20)
%            tol   -- tolerance for Newton iteration: norm of gradient (1e-6)
%            cgits -- max. number of BiCGstab iterations (1000)
%            cgtol -- tolerance for BiCGstab (1e-6)
%            verb  -- verbosity (0: no display)
% Output: u      -- reconstruction
%
% Christian Clason (christian.clason@uni-graz.at)
% Bangti Jin       (btjin@math.tamu.edu)
% February 27, 2011

%% Setup parameters
nel   = size(fem.xc,1);     % problem size: number of elements

% Default parameters
if nargin <= 5 || isempty(opts)
    maxit = 20;             % max. number of semi-smooth Newton iterations
    tol   = 1e-6;           % tolerance for Newton iteration (norm of gradient)
    cgits = 1000;           % max. number of BiCGstab iterations
    cgtol = 1e-6;           % tolerance for BiCGstab
    verb  = 1;              % verbosity (0,1: no display, 2: only path-following, 3: also SSN)
else
    maxit = opts.maxit;
    tol   = opts.tol;
    cgits = opts.cgits;
    cgtol = opts.cgtol;
    verb  = opts.verb;
end
% creating Mf as right hand side

P   = fem.P;      % Projection onto piecewise constants for linear functions
M   = fem.M;      % Mass matrix for piecewise constant on boundary
%Y   = fem.Y;      % Y(v)*du assembles weighted mass matrix for piecewise constant du and linear v
A   = fem.A;                 % differential operator
obs = fem.obs;    % restriction to nodes on observation boundary

%% SSN loop

% initialize iterates
u = ones(nel,1);
iter = 0;

while (iter<maxit)
    
    %catch notspd %#ok<NASGU>
        S = @(f) A\f;       % fallback if numerically semidefinite
    %end
    
    y = S(M*Pu(fem,u)+Mf);                  % state
    p = 0*y;
    p(obs) = y(obs)-yd(obs);                  % residual
%    
    q = S(-(M*p));              % adjoint
%    
%    My = Y(y);
%    Mq = Y(q);         % precompute mass matrices for Hessian application
%   
    % gradient
    rhs  = -(alpha*u + P(y.*q));
    nres = norm(rhs)/sqrt(nel);    
    % Hessian
    C = @(x) applyHess(x,S,P,M,obs,alpha);
    
    % Newton step
    [dx,flag] = bicgstab(C,rhs,cgtol,cgits); %#ok<NASGU>
    
    u = u + dx;
    
    if verb > 1
        fprintf('Newton iteration %d: |grad|=%e, |u-uex|=%e \n',iter,nres); % norm(u-ue)/sqrt(nel)
    end
    
    if nres<tol    % terminate if active sets coincide and norm of gradient small
        break
    end
    
    iter   = iter+1;
    
end
%end main function
end
%% Application of Hessian
function Hdu = applyHess(x,S,P,M,obs,alpha)

dy = S(-M*x);                 % differential change in state
dp = 0*dy; dp(obs) = dy(obs);
dq = S(M*dp);            % differential change in adjoint

Hdu = alpha*x + P(dq);
%end function applyHess
% function for projecting u
end
%%
function uu=Pu(fem,u)
    uuu=[0;u];
    uu=uuu;
    for i=1:size(fem.xx,1)-1
        uu =[uu;uuu]; 
    end 
end
