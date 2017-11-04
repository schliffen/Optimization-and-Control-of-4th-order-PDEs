% -------------
% constraint 4th order optimal control problem
% ------------
clear all
close all
n = input(' Enter the size of the problem~ n :   ')
%
%
fem = setFEM2d(n);
% setting parameters
ny = size(fem.L,2);
nu = size(fem.M,2);
tol=0.0001;
c  = 1;
vu = .4;
itmax = 4;
% Loding the proper matrices (this is for boundary control problem)
  M = fem.L;  % mass matrix in domain
  Bu= fem.M;  % mass matrix on boundary
  By= fem.M;  % mass matrix on boundary 
  L = fem.A;  % stiffness matrix with boundary condition
% constant parameters
% constant vectors
kk=size(fem.xx,2); 
yd = 0.2*ones(ny,1);  yd(1:16*kk,1) = 0.002;yd(end-16*kk+1:end,1) = 0.002;
%yd = 0*ones(ny,1);  yd(1:ny/2,1) = 0.02;%yd(end-16*kk+1:end,1) = 0.002;
% creation of target function

%    yd=zeros(kk,1);
%     for k=1:size(fem.xx,2)
%       for kk=1:size(fem.yy,2)
%   if (fem.xx(1,k)-.5)^2+ (fem.yy(1,kk)-.5)^2 <=.3
%       yd(kk,k) =5;
%   else
%       yd(kk,k)=0;
%   end
%       end
%     end
%     yd = reshape(yd,ny,1);
  % Control bounds
  ua = .00001*ones(nu,1);   % untere Grenze (punktweise)
  ub = 10*ones(nu,1); % obere Grenze (punktweise)
%  
for ra=1:5
%
  y = zeros(ny,1); u = 20*ones(nu,1);  % same # of dofs (distributed control)
  p = zeros(ny,1); mu = zeros(nu,1);
%
% Initialize flags and counters
iter = 0; done = 0;
stop = false;
Y = [y;u;p;mu];
% for comparing active sets
Apl = zeros(nu,1) ; Amin = zeros(0,1);
% Do a simple semismooth Newton loop
   bqr=zeros(ny,1); bqr(ny/2-50:ny/2+50,1)=.01;
   %bqr(1:nu/2,1)=.001*rand(nu/2,1);
   bqr=bqr+0.001*randn(ny,1);
   while (~stop)
%         % Determine active and inactive sets
  % ---------------------------------------------------
  Aplus  = find(mu + c * (u - ub) > 0);
  Aminus = find(mu + c * (u - ua) < 0);
  
  % Solve the Newton system by Bramble-Pasciak like cg
  % ---------------------------------------------------
% computing objective function   
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
%  P = fem.P;  % projector function
  % Set up the saddle point blocks
  A = [M  Z1; Z1' vu*Bu];
  B = [L  -Bu; Z2' PA];
  C = [L' Z2; -By' PA'];
  %
     % Set up right hand side
  bx = [M*yd; zeros(nu,1)];
  bq = zeros(nu,1);
   if size(Aplus,1)>1
       bq(Aplus) = ub(Aplus);
   end 
  bq(Aminus) = ua(Aminus);
  bq = [bqr; PA*bq];
  %
% constructing matrices
  K = [A B';B zeros(ny+nA,ny+nA)]; 
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
%  dY = K\b;
%
%P=[A zeros(size(A));zeros(size(A)) B];
[dY,flag,relres,cgiter] = bicgstabl(K,b,tol,30);
  time_fmincon = toc(startTime);
clear K b  
%Y(1:2*ny+nu) = Y(1:2*ny+nu)+dY(1:2*ny+nu);mu(union(Aplus,Aminus)) +
%u =  dY(ny+1:ny+nu); %
  mu(union(Aplus,Aminus)) = mu(union(Aplus,Aminus)) + dY(2*ny+nu+1:end);
%          
% iteration updating
           iter = iter + 1;   
% Checking stopping criteria
 ntol=norm(dY(1:ny+nu),2)
  		if (iter >= itmax)||ntol<tol%size(Apl,1) == 0 && size(Amin,1) == 0 && size(Aminus,1)==0 && size(Aplus,1)== 0  
			stop = true;
        end

%saving last control 
%update iteration
if ~stop
    if iter == 1
    Ly=dY;
    else
    Ly(1+ny:2*ny+nu,1) = Ly(1+ny:2*ny+nu,1)+ dY(1+ny:2*ny+nu,1);  
    end
elseif iter == 1
     Ly=dY;
end
u = Ly(ny+1:ny+nu); 
%mu(union(Aplus,Aminus)) = Ly(2*ny+nu+1:end);
% Clearing active sets   
%Apl = Aplus; Amin = Aminus;
           clear Aplus Aminus;
%y=dY(1:ny);
%f=.5*(y-yd)'*M*(y-yd)+.5*vu*u'*Bu*u 
end
if ra>1
    Lr=[Lr Ly(1:ny+nu)];
else
    Lr=Ly(1:ny+nu);
end
clear Ly u y mu dY p mu
end
y=Lr(1:ny,:);
%u=Ly(ny+1:ny+nu);
%p=Ly(ny+nu+1:2*ny+nu);
fprintf('back slash solution takes %g seconds.\n',...
time_fmincon); % Ausgabe der benoetigten Zeit
% Postprocessing
%objective function
%u=P(u);
x = linspace(0,1,n)';
%figure(1)
u=sum(Lr(ny+1:ny+nu,:),2)/ra;
y=sum(Lr(1:ny,:),2)/ra;
uu=reshape(u,n,n);
%yy=reshape(y,n,n);
%plot(x,uu(:,1),'b')
%hold on
%uu=mean(uu,2);
%trisurf(fem.tri,fem.xx,fem.yy,0*xx,y,'edgecolor','k','facecolor','interp');
% proper ploting
%xcord=fem.xx(2:end-1,:);
%ycord=fem.yy(2:end-1,:);
%finy=y(n:ny-n,1)
%
figure(1);plot(x,uu(:,1),'r');
figure(5);plot(x,uu(:,n),'r');
figure(2);trisurf(fem.tri,fem.xx,fem.yy,y,'edgecolor','k','facecolor','interp')
figure(3);trisurf(fem.tri,fem.xx,fem.yy,bqr,'edgecolor','k','facecolor','interp')
figure(4);trisurf(fem.tri,fem.xx,fem.yy,M*yd,'edgecolor','b','facecolor','interp')
f=.5*(y-yd)'*M*(y-yd)+.05*vu*u'*Bu*u
var=zeros(ny,1);
for i=1:ra
    var=var + (Lr(1:ny,i)-y).^2;
end
figure(6);trisurf(fem.tri,fem.xx,fem.yy,var,'edgecolor','k','facecolor','interp')
E = abs(L*y-Bu*u-bqr);
figure(7);trisurf(fem.tri,fem.xx,fem.yy,.0001*E,'edgecolor','k','facecolor','interp')
%tracking error
Et=abs(y-yd);
figure(8);trisurf(fem.tri,fem.xx,fem.yy,Et,'edgecolor','k','facecolor','interp')