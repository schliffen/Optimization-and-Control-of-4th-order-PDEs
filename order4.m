% -------------
% constraint 4th order optimal control problem
% ------------
n = input('enter the size of the problem~ n = ')
%mat = material;
%X = netz(n); % Diskretisierung erstellen (Stuetzwerte 1SP, Steuerung 2SP)
ua = -ones(n,1); % untere Grenze (punktweise)
ub = ones(n,1); % obere Grenze (punktweise)
%
%[L,M,Bu,By]=statO4(n);
%
fem = setFEM2d(n);
% setting parameters
ny = size(fem.Y,2);
nu = size(fem.M,2);
c = 1;
vu = .1;
tol=.001;
yd = ones(ny,1);
%
y = zeros(ny,1); u = ones(nu,1);  % same # of dofs (distributed control)
p = zeros(ny,1); mu = zeros(nu,1);
%
% Initialize flags and counters
iter = 0; done = 0;
stop = false;
itmax = 100;
Y = [y;u;p;mu];
% Do a simple semismooth Newton loop
while (~stop)
%    
     % Determine active and inactive sets
  % ---------------------------------------------------
  Aplus  = find(mu + c * (u - ub) > 0);
  Aminus = find(mu + c * (u - ua) < 0);
  
  % Solve the Newton system by Bramble-Pasciak like cg
  % ---------------------------------------------------
  % Set up the active set projector
  PA = zeros(nu,1);
  PA(Aplus) = 1; PA(Aminus) = 1;
  PA = spdiags(PA,0,nu,nu);
  PA = PA(union(Aplus,Aminus),:);
  nA = length(Aplus)+length(Aminus);
  fprintf(' %4d  %6d  %6d',iter,length(Aplus),length(Aminus));
  %
  % Set up some zero matrices
  Z1 = sparse(ny,nu); Z2 = sparse(ny,nA);
  %
  % Loding the proper matrices (this is for boundary control problem)
  M = fem.Y;
  Bu= fem.M;
  By= fem.M;
  L = fem.A;
  % Set up the saddle point blocks
  A = [M  Z1; Z1' vu*Bu];
  B = [L  -By; Z2' PA];
  C = [L' Z2; -By' PA'];
  %
     % Set up right hand side
  bx = [M*yd; zeros(nu,1)];
  bq = zeros(nu,1);
  if size(Aplus,1)>1
  bq(Aplus) = ub(Aplus);
  end
  bq(Aminus) = ua(Aminus);
  bq = [zeros(ny,1); PA*bq];
% constructing matrices
K = [A C;B zeros(ny+nA,ny+nA)]; 
b = [bx;bq];
% Optimierungsparameter festlegen
%options=optimset(@fmincon); % urspruengliche Parameter
startTime=tic; % Zeitmessung fuer die Optimierung
%options=optimset(options,'Display','iter','Algorithm','active-set',...
%'TolFun',1e-10,'TolCon',1e-10,'TolX',1e-10,...
%'UseParallel','always','GradObj','on',...
%'DerivativeCheck','off','Diagnostics','on',...
%'FunValCheck','on')
% Optimierung starten
dY = K\b;
%
time_fmincon = toc(startTime);
%Y(1:2*ny+nu) = Y(1:2*ny+nu)+dY(1:2*ny+nu);
u = dY(ny+1:ny+nu); 
mu(union(Aplus,Aminus)) = mu(union(Aplus,Aminus)) + dY(2*ny+nu+1:end);
iter = iter + 1;
condd=norm(dY(1:ny)-y);
		if (iter >= itmax)|| (condd<tol)
			stop = true;
        end  
y=dY(1:ny);
f=.5*(y-yd)'*M*(y-yd)+.5*vu*u'*Bu*u        
end
y=dY(1:ny);
u=dY(ny+1:ny+nu);
p=dY(ny+nu+1:2*ny+nu);
fprintf('back slash solution takes %g seconds.\n',...
time_fmincon); % Ausgabe der benoetigten Zeit
% Postprocessing
%objective function 
f=.5*(y-yd)'*M*(y-yd)+.05*u'*Bu*u

