% ut + vel*ux = 0 
% Dirichlet u(0,t) = 0
% or periodic BC 
% Solved with DG weak form (arbitrary poly. order) and various fluxes
function A1Driver_1D_waveEq
clear all
clc

% Domain bounds
x_min=0.0
x_max=2*pi

% Advection velocity
vel=1.0

% Check errors with periodic BC
% try NumElements=10 and N = 2 --> DoF=(P+1)*Nel=30  --> Error=0.1312 (for 1 iteration)
% try NumElements=10 and N = 4 --> DoF=(P+1)*Nel=50  --> Error=0.1005
% try NumElements=10 and N = 6 --> DoF=(P+1)*Nel=70  --> Error=0.084854

% Number of elements in the mesh 
NumElements=10 % min number of elements is 2

% Order of polymomials used for approximation 
N = 6

% scheme  = numerical flux
%         = 0 (Central-Diff) --> Unstable
%         = 1 (Lax-Friedrichs)
%         = 2 (Lax-Wendroff)
%         = 3 (Murman-Roe)
%         = 4 (Godunov)
%         = 5 (Roe with entropy fix)
scheme=1;

%% Alternative flux formulation
% Upwind or central fluxes
% 0.0 --> Upwind
% 1.0 --> Central
%alfaFlux=0.0; 

% these parameters control what type of Jabobi polynomials: alpha=0.0; beta=0.0 is Legendre
alpha=0.0; beta=0.0; 

% Euler explicit
nt=80  % total time steps
dt=0.01 % time step

delay = 0.1; % frame delay for movie

% Boundary Condition: 
% 1 --> Periodic
% 0 --> Dirichlet
BC=1


%%%%%%%%%%%%%%%
%MESH
%%%%%%%%%%%%%%%
deltaX=(x_max-x_min)/NumElements; 


% define X coord in eahc mesh node
[x,w] = JacobiGQ(alpha,beta,N);   % if the polynomial is always the same this can be taken outside

% define node structure
for i=1:NumElements  
         
    mesh(i).x_left=x_min+(i-1)*deltaX;
    mesh(i).x_right=x_min+(i)*deltaX;
        
    mesh(i).xglobal=zeros(N+1,1);
    mesh(i).xglobal=Mapp_r2x(mesh(i).x_left,mesh(i).x_right,x); % this is a vector with all Gauss poijt for each element    
  
    mesh(i).normal_left=1.0;  % -o--> el <--o-
    mesh(i).normal_right=-1.0;
    
    mesh(i).flux_left=zeros(N+1,1);  
    mesh(i).flux_right=zeros(N+1,1);  
    
    
end

% constant Jacobian (same for all element here), but can be changed here
J=ComputeJacobian(mesh(1).x_left,mesh(1).x_right);

for i=1:NumElements  
         
    mesh(i).u=zeros(N+1,1);
    mesh(i).u=Project_phys2modal(alpha,beta,N,J,sin(2.0*mesh(i).xglobal)); % this is initial condition
    
%    this fields could be included in each elements
%    mesh(i).local_poly=N;
%    mesh(i).weight=w;
%    mesh(i).f=zeros(N+1,1); 
%    mesh(i).J=0.0; 
end

% this is used to identify nodes when computing fluxes at nodes 
Left_node=-1;
Right_node=1; 
%%%%%%%%%%%%%%%
% Construct local to element matrices
%%%%%%%%%%%%%%%
%M = GenerateMass_matrix1D(alpha,beta,N); % this one should be the identity (sanity check)
D = GenerateDeriv1D(alpha,beta,N);
%DD = Generate2ndDeriv1D(alpha,beta,N);

%Scale matrices with Jacobian of the element (here all the same because uniform poly and equidistant nodes)
%M=M*J; % not used because M=Id
%DD=DD/J;

%%%%%%%%%%%%%%%
% Construct global matrix
%%%%%%%%%%%%%%%
dim_global=NumElements*(N+1); % overall problem dimension is num element x number modes per element

% allocate global matrices
Uglobal=zeros(dim_global,1);
Xplot=zeros(dim_global,1);
Uplot=zeros(dim_global,1);
Fluxglobal=zeros(dim_global,1);  
Matrixglobal=zeros(dim_global,dim_global);


% add local matrices to the diagonal of the golbal matrix
for i=1:NumElements  
   loc_start=1+(i-1)*(N+1);
   loc_end=(i)*(N+1);
   Matrixglobal(loc_start:loc_end,loc_start:loc_end)=vel*D(:,:);
end

% spy(Matrixglobal)

%%%%%%%%%%%%%%%
% Initial plot (will be updated later)
%%%%%%%%%%%%%%%

OnlyNodes=zeros(NumElements+1,1);
OnlyUatNodes=zeros(NumElements+1,1);
%
count=1;
for i=1:NumElements
    OnlyNodes(i)=mesh(i).x_left;  % only h-mesh nodal information
    OnlyUatNodes(i)=Interpol_modal2phy(alpha,beta,N,mesh(end).u,Left_node);
    for k=1:N+1  
       Xplot(count)=mesh(i).xglobal(k) % info at all h/p nodes (including Gauss points)
       
       Uglobal(count)=mesh(i).u(k);
       
       Uplot(count)=Interpol_modal2phy(alpha,beta,N,mesh(i).u,x(k))
      
       count=count+1;
    end
end
OnlyNodes(NumElements+1)=mesh(NumElements).x_right;
OnlyUatNodes(NumElements+1)=Interpol_modal2phy(alpha,beta,N,mesh(NumElements).u,Right_node);

plot(OnlyNodes,OnlyUatNodes,'bo'); hold('on'); 
for i=1:NumElements  
   plot(mesh(i).xglobal , Interpol_modal2phy_all_points(alpha,beta,N,mesh(i).u),'g-'); 
end
%  data will be updated later
p = plot(Xplot,Uplot,'r');
legend('Nodal values', 'Initial Condition','DG solution');
axis([x_min x_max -1.5 1.5]); hold('off'); 
%you can use polyval to plot a polynomial and not lines between high order nodes, but movie is difficult to make

%%%%%%%%%%%%%%%
% Time marching (using Euler explicit)
%%%%%%%%%%%%%%%

CFL=vel*dt/(Mapp_r2x(mesh(1).x_left,mesh(1).x_right,x(2))-Mapp_r2x(mesh(1).x_left,mesh(1).x_right,x(1)))

dx=vel*dt/CFL;

for k=1:nt  
display('iteration: '); k

%%%%%%%%%%%%%%%
% Update flux vector, including boundary conditions
%%%%%%%%%%%%%%%
for i=2:NumElements-1  % this hould be nodes!! not elements 
  
   % LEFT NODE OF ELEMENT i where  fnode_ext--o--fnode_int        
    fnode_ext=Interpol_modal2phy(alpha,beta,N,vel*mesh(i-1).u,Right_node); % evaluate flux at element edge (-1 or 1)      
    fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(i).u,Left_node);            

%   mesh(i).flux_left=ComputeFlux_Lax1D(alpha,beta,N,fnode_ext,fnode_int,mesh(i-1).normal_right,mesh(i).normal_left,Left_node,vel,alfaFlux); % compute flux at left (x=-1) node of element i
    mesh(i).flux_left=ComputeFluxes(alpha,beta,N,fnode_ext,fnode_int,Left_node,vel,scheme,dt,dx); % compute flux at left (x=-1) node of element i 
  
   % RIGTH NODE OF ELEMENT i where:  fnode_int--o--fnode_ext          
    fnode_ext=Interpol_modal2phy(alpha,beta,N,vel*mesh(i+1).u,Left_node);  
    fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(i).u,Right_node); % evaluate flux at element edge (-1 or 1)             

%   mesh(i).flux_right=ComputeFlux_Lax1D(alpha,beta,N,fnode_int,fnode_ext,mesh(i).normal_right,mesh(i+1).normal_left,Right_node,vel,alfaFlux); % compute flux at right (x=+1) node of element i
    mesh(i).flux_right=ComputeFluxes(alpha,beta,N,fnode_int,fnode_ext,Right_node,vel,scheme,dt,dx); % compute flux at right (x=+1) node of element i
end

% Dirichlet
if(BC==0)
if (k==1) 
    display('Dirichlet BC at x_min') 
end

% FIRST ELEMENT
%BC LEFT NODE  OF ELEMENT 1: x=0
    %u_ext=Project_phys2modal(alpha,beta,N,J,-sin(vel*k*dt)*ones(N+1,1));%zeros(N+1,1);
    u_ext=Project_phys2modal(alpha,beta,N,J,zeros(N+1,1)); % ghost node
    fnode_ext=Interpol_modal2phy(alpha,beta,N,u_ext,Right_node); % ghost element with zeros  
    normal_temp=1.0;
    fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(1).u,Left_node);  
 
%   mesh(1).flux_left=ComputeFlux_Lax1D(alpha,beta,N,fnode_ext,fnode_int,normal_temp,mesh(1).normal_left,Left_node,vel,alfaFlux); 
    mesh(1).flux_left=ComputeFluxes(alpha,beta,N,fnode_ext,fnode_int,Left_node,vel,scheme,dt,dx); 

% RIGTH NODE OF ELEMENT 1         
    fnode_ext=Interpol_modal2phy(alpha,beta,N,vel*mesh(2).u,Left_node); 
    fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(1).u,Right_node); % evaluate flux at element edge (-1 or 1)            
 
%   mesh(1).flux_right=ComputeFlux_Lax1D(alpha,beta,N,fnode_int,fnode_ext,mesh(1).normal_right,mesh(2).normal_left,Right_node,vel,alfaFlux); % compute flux at right (x=+1) node of element i
    mesh(1).flux_right=ComputeFluxes(alpha,beta,N,fnode_int,fnode_ext,Right_node,vel,scheme,dt,dx); % compute flux at right (x=+1) node of element i
 
% LAST ELEMENT  
% LEFT NODE OF ELEMENT NumElements        
    fnode_ext=Interpol_modal2phy(alpha,beta,N,vel*mesh(NumElements-1).u,Right_node); % evaluate flux at element edge (-1 or 1)      
    fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(NumElements).u,Left_node);            

%   mesh(NumElements).flux_left=ComputeFlux_Lax1D(alpha,beta,N,fnode_ext,fnode_int,mesh(NumElements-1).normal_right,mesh(NumElements).normal_left,Left_node,vel,alfaFlux); % compute flux at left (x=-1) node of element i
    mesh(NumElements).flux_left=ComputeFluxes(alpha,beta,N,fnode_ext,fnode_int,Left_node,vel,scheme,dt,dx); % compute flux at left (x=-1) node of element i
 
% RIGHT NODE OF ELEMENT NumElements     
    fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(NumElements).u,Right_node); % evaluate flux at element edge (-1 or 1)          
    %u_ext=zeros(N+1,1);
    %fnode_ext=Interpol_modal2phy(alpha,beta,N,u_ext,Left_node); % ghost element with zeros  
    fnode_ext=fnode_int;
    normal_temp=1.0; 
 
%   mesh(NumElements).flux_right=ComputeFlux_Lax1D(alpha,beta,N,fnode_int,fnode_ext,mesh(NumElements).normal_right,normal_temp,Right_node,vel,alfaFlux); % compute flux at right (x=+1) node of element i
    mesh(NumElements).flux_right=ComputeFluxes(alpha,beta,N,fnode_int,fnode_ext,Right_node,vel,scheme,dt,dx); % compute flux at right (x=+1) node of element i

end

% Periodic
if(BC==1)  
if (k==1) 
    display('Periodic BC') 
end

%FIRST ELEMENT:
% LEFT NODE OF ELEMENT 1
    fnode_ext=Interpol_modal2phy(alpha,beta,N,vel*mesh(NumElements).u,Right_node); 
    fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(1).u,Left_node);  

%   mesh(1).flux_left=ComputeFlux_Lax1D(alpha,beta,N,fnode_ext,fnode_int,mesh(NumElements).normal_right,mesh(1).normal_left,Left_node,vel,alfaFlux);
    mesh(1).flux_left=ComputeFluxes(alpha,beta,N,fnode_ext,fnode_int,Left_node,vel,scheme,dt,dx);
  
% RIGTH NODE OF ELEMENT 1             
    fnode_ext=Interpol_modal2phy(alpha,beta,N,vel*mesh(2).u,Left_node); 
    fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(1).u,Right_node);            

%   mesh(1).flux_right=ComputeFlux_Lax1D(alpha,beta,N,fnode_int,fnode_ext,mesh(1).normal_right,mesh(2).normal_left,Right_node,vel,alfaFlux); % compute flux at right (x=+1) node of element i
    mesh(1).flux_right=ComputeFluxes(alpha,beta,N,fnode_int,fnode_ext,Right_node,vel,scheme,dt,dx); % compute flux at right (x=+1) node of element i

%LAST ELEMENT: 
% LEFT NODE OF ELEMENT NumElements        
    fnode_ext=Interpol_modal2phy(alpha,beta,N,vel*mesh(NumElements-1).u,Right_node); % evaluate flux at element edge (-1 or 1)      
    fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(NumElements).u,Left_node);            
 
%   mesh(NumElements).flux_left=ComputeFlux_Lax1D(alpha,beta,N,fnode_ext,fnode_int,mesh(NumElements-1).normal_right,mesh(NumElements).normal_left,Left_node,vel,alfaFlux); % compute flux at left (x=-1) node of element i
    mesh(NumElements).flux_left=ComputeFluxes(alpha,beta,N,fnode_ext,fnode_int,Left_node,vel,scheme,dt,dx); % compute flux at left (x=-1) node of element i
   
% RIGHT NODE OF ELEMENT NumElements         
   fnode_ext=Interpol_modal2phy(alpha,beta,N,vel*mesh(1).u,Left_node); % evaluate flux at element edge (-1 or 1)      
   fnode_int=Interpol_modal2phy(alpha,beta,N,vel*mesh(NumElements).u,Right_node);            

%   mesh(NumElements).flux_right=ComputeFlux_Lax1D(alpha,beta,N,fnode_int,fnode_ext,mesh(NumElements).normal_right,mesh(1).normal_left,Right_node,vel,alfaFlux); % compute flux at left (x=-1) node of element i
    mesh(NumElements).flux_right=ComputeFluxes(alpha,beta,N,fnode_int,fnode_ext,Right_node,vel,scheme,dt,dx); % compute flux at left (x=-1) node of element i
end %(periodicBC==1)   
 
 
% Update Flux vector and vector to plot solution
count=1;
for i=1:NumElements  
  for k=1:N+1
     Fluxglobal(count)=-mesh(i).flux_right(k)+mesh(i).flux_left(k); % update fluxes in golbal vector fluxes
     count=count+1;
  end
end 
%%%%%%%%%%%%%%

Uglobal=Uglobal+dt/J*(Matrixglobal'*Uglobal+Fluxglobal) ; % 1/J comes from the mass mattrix

count=1;
for i=1:NumElements  
  for k=1:N+1  
      

      mesh(i).u(k)=Uglobal(count)  % copy new data to nodal strucuture 
    
      Uplot(count)=Interpol_modal2phy(alpha,beta,N,mesh(i).u,x(k))
     count=count+1;
  end
end  


% redraw solution
set(p,'ydata',Uplot);
pause(delay)
 
end % time marching


count=1;
for i=1:NumElements  
  for k=1:N+1  
      
      Uexact(count)=sin(2.0*(mesh(i).xglobal(k)-vel*dt*nt));
      count=count+1;
  end
end 

hold on; plot(Xplot,Uexact,'b'); hold off;
%legend('Nodal values', 'Initial Condition','DG solution', 'Exact solution')

error=(1/dim_global)*norm(Uexact(:,1)-Uplot(:,1))%/norm(Uexact(:,1))
%error=trapz(abs(Uexact(:,1)-Uplot(:,1))) 




count=1;
for i=1:NumElements  
  for k=1:N+1  
      Nodos(count) = mesh(i).xglobal(k)
      count=count+1;
  end
end 

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSTRUCCIÓN DE LOS ARCHIVOS DE TEXTO

Nodos_datos=[0.0160    0.0812    0.1867    0.3142    0.4417    0.5471    0.6123    0.6443    0.7095    0.8150    0.9425    1.0700    1.1754    1.2406    1.2726    1.3378    1.4433    1.5708   1.6983    1.8038    1.8690    1.9009    1.9662    2.0716    2.1991    2.3266    2.4321    2.4973    2.5293    2.5945    2.6999    2.8274    2.9549    3.0604    3.1256    3.1576    3.2228    3.3283    3.4558    3.5833    3.6887    3.7539    3.7859    3.8511    3.9566    4.0841    4.2116    4.3170    4.3822    4.4142    4.4794    4.5849    4.7124    4.8399     4.9453    5.0106    5.0425    5.1077    5.2132    5.3407    5.4682    5.5737    5.6389    5.6709    5.7361    5.8415    5.9690    6.0965    6.2020    6.2672]
fileID = fopen('Nodos.txt','w');
fprintf(fileID,'%6s %12s\n','Nodos');
fprintf(fileID,'%6.2f %12.8f\n',Nodos_datos);
fclose(fileID);

count=1
for i=0:69
    tiempo(count) = i +1
    count=count+1
end


U_datos = [Uplot(:,1), Nodos(:), tiempo(:) ]
fileID = fopen('U.txt','w');
fprintf(fileID,'%6s %12s\n','U');
fprintf(fileID,'%6.2f %12.8f\n',U_datos);
fclose(fileID);

end % end main function












%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% Additional functions 
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

function [P] = JacobiP(x,alpha,beta,N)
% function [P] = JacobiP(x,alpha,beta,N)
% Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
%          (alpha+beta <> -1) at points x for order N and returns P[1:length(xp))]
% Note   : They are normalized to be orthonormal.

% Turn points into row if needed.
xp = x; dims = size(xp);
if (dims(2)==1) xp = xp';  end;

PL = zeros(N+1,length(xp)); 

% Initial values P_0(x) and P_1(x)
gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
    gamma(beta+1)/gamma(alpha+beta+1);
PL(1,:) = 1.0/sqrt(gamma0);
if (N==0) P=PL'; return; end;
gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
PL(2,:) = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1);
if (N==1) P=PL(N+1,:)'; return; end;

% Repeat value in recurrence.
aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

% Forward recurrence using the symmetry of the recurrence.
for i=1:N-1
  h1 = 2*i+alpha+beta;
  anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*...
      (i+1+beta)/(h1+1)/(h1+3));
  bnew = - (alpha^2-beta^2)/h1/(h1+2);
  PL(i+2,:) = 1/anew*( -aold*PL(i,:) + (xp-bnew).*PL(i+1,:));
  aold =anew;
end;

P = PL(N+1,:)';
end

%%%%%%%%%%%%%%%

function [x,w] = JacobiGQ(alpha,beta,N)
% function [x,w] = JacobiGQ(alpha,beta,N)
% Purpose: Compute the N'th order Gauss quadrature points, x, 
%          and weights, w, associated with the Jacobi 
%          polynomial, of type (alpha,beta) > -1 ( <> -0.5).

if (N==0) x(1)= -(alpha-beta)/(alpha+beta+2); w(1) = 2; return; end;

% Form symmetric matrix from recurrence.
J = zeros(N+1);
h1 = 2*(0:N)+alpha+beta;
J = diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1) + ...
    diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta).*...
    ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);
if (alpha+beta<10*eps) J(1,1)=0.0;end;
J = J + J';

%Compute quadrature by eigenvalue solve
[V,D] = eig(J); x = diag(D);
w = (V(1,:)').^2*2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
    gamma(beta+1)/gamma(alpha+beta+1);
end

%%%%%%%%%%%%%%%

function [dP] = GradJacobiP(r, alpha, beta, N)
% function [dP] = GradJacobiP(r, alpha, beta, N);
% Purpose: Evaluate the derivative of the Jacobi polynomial of type (alpha,beta)>-1,
%	       at points r for order N and returns dP[1:length(r))]        

dP = zeros(length(r), 1);
  if(N == 0)
    dP(:,:) = 0.0; 
  else
    dP = sqrt(N*(N+alpha+beta+1))*JacobiP(r(:),alpha+1,beta+1, N-1);
  end;
end


%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% Generate Mass matrix (here we have the M=Id, this is used for the checking implemetation)
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

function [M] = GenerateMass_matrix1D(alpha,beta,N)
% function [D] = Mass amtrix
% integral of p*p dr

[x,w] = JacobiGQ(alpha,beta,N);
M = zeros((N+1),(N+1));

% Initialize matrix
for i=1:N+1
  P1 = JacobiP(x, alpha, beta, i-1); 
  for j=1:N+1
   P2 = JacobiP(x,alpha,beta,j-1);
       for k=1:N+1
         M(i,j) =M(i,j)+w(k)*P1(k)*P2(k) ; 
       end
   end
end

end

%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% Generate stiffness matrix (first order derivative)
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

function [D] = GenerateDeriv1D(alpha,beta,N)
% function [D] = First derivative matrix
% integral of dp/dr*p dr

[x,w] = JacobiGQ(alpha,beta,N);
D = zeros((N+1),(N+1));

% Initialize matrix
for i=1:N+1
   P = JacobiP(x,alpha,beta,i-1); 
  for j=1:N+1
       dP = GradJacobiP(x, alpha, beta, j-1);
     for k=1:N+1
       D(i,j) =D(i,j)+w(k)*P(k)*dP(k);
     end
   end
end
end 

%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% Generate Laplacian matrix (second order derivative)
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

function [DD] = Generate2ndDeriv1D(alpha,beta,N)
% function [DD] = Second derivative matrix
% integral of dp/dr*dp/dr dr

[x,w] = JacobiGQ(alpha,beta,N);
DD = zeros((N+1),(N+1));

% Initialize matrix
for i=1:N+1
  dP1 = GradJacobiP(x, alpha, beta, i-1); 
  for j=1:N+1
    dP2 = GradJacobiP(x,alpha,beta,j-1);
     for k=1:N+1
        DD(i,j) =DD(i,j)+w(k)*dP1(k)*dP2(k) ;
     end
  end
end
end

%%%%%%%%%%%%%%%
% Project physical value to modal (i.e. calculate modal coeffs.)
%%%%%%%%%%%%%%%

function modal=Project_phys2modal(alpha,beta,N,J,f)
  % the operation is: u_modal=M^-1 \int u_phy*phi --> M^1 is 1/J here
  % here the mass matrix is the identity (if not this needs to be change)

[x,w] = JacobiGQ(alpha,beta,N);
modal=zeros((N+1),1) ;

  for i=1:N+1  
   P = JacobiP(x, alpha, beta, i-1); 
     for k=1:N+1
      modal(i,1) =modal(i,1)+w(k)*P(k)*f(k); %*J/J ; 
     end
  end   
end 

%%%%%%%%%%%%%%%
% Interpolate modal values to physical values, for one particular point (e.g. boundary values)
%%%%%%%%%%%%%%%

function Phy=Interpol_modal2phy(alpha,beta,N,f_modal,at_x)
  % the operation is: u_phy=Sum(u_modal*phi)
Phy=0.0;

  for i=1:N+1  
   P = JacobiP(at_x, alpha, beta, i-1); 
      Phy =Phy+f_modal(i)*P ; %modes start at 1 and poly array at 1, but min order is zero
  end
end 

%%%%%%%%%%%%%%%
% Interpolate modal values to physical values, for all element Gauss points
%%%%%%%%%%%%%%%

function Phy=Interpol_modal2phy_all_points(alpha,beta,N,f_modal)
  % the operation is: u_phy=Sum(u_modal*phi)
Phy=0.0;
[x,w] = JacobiGQ(alpha,beta,N);

  for i=1:N+1  
   P = JacobiP(x, alpha, beta, i-1); 
      Phy =Phy+f_modal(i)*P ; %modes start at 1 and poly array at 1, but min order is zero
  end
end 

%%%%%%%%%%%%%%%
% Map physical to computational space
%%%%%%%%%%%%%%%

function xglobal=Mapp_r2x(a,b,xlocal);
 % 1D affine mappings r=[-1,1]-->x=[x_coord_left,x_coord_right]
  h=abs(a-b);
  xglobal=a+(1.0+xlocal)*(h)/2.0; % r runs from -1 to 1
end

function xlocal=UnMapp_x2r(a,b,xglobal);
 % 1D affine mappings r=[-1,1]<--x=[x_coord_left,x_coord_right]
  h=abs(a-b);
  xlocal=(xglobal-a)*2.0/h-1.0; % r runs from -1 to 1
end

%%%%%%%%%%%%%%%
% Compute element Jacobian dr/dx
%%%%%%%%%%%%%%%

function J=ComputeJacobian(a,b);
  J=abs(a-b)/2.0;
end


%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% Alternative Flux
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

function Flux=ComputeFlux_Lax1D(alpha,beta,N,f1,f2,normal2,normal1,at_x,vel,alfaFlux)

Flux=zeros((N+1),1);

temp=((f1+f2)/2.0+abs(vel)*(f1*normal1+f2*normal2)*(1.0-alfaFlux)/2.0);
  for i=1:N+1  
   P = JacobiP(at_x, alpha, beta, i-1);   
      Flux(i,1) =temp*P ;      
  end   
end



%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%------------------------------------------------------------------------------
% Flux of conservation law
%------------------------------------------------------------------------------
function flux = FLUX(u)
%flux = 0.5*u^2; %Burgers
flux = u; % linear advection
end

%------------------------------------------------------------------------------
% Derivative of flux in conservation law
%------------------------------------------------------------------------------
function dflux = DFLUX(u)
%dflux = u;  % Burger's
dflux = 1.0; % linear advection
end

%------------------------------------------------------------------------------
% Approximate Derivative of flux in conservation law
%------------------------------------------------------------------------------
function a = ADFLUX(u,v)
if abs(u-v) > 1.e-14
   a = (FLUX(u) - FLUX(v))/(u-v);
else
   a = DFLUX(u);
end
end

%------------------------------------------------------------------------------
% Numerical flux function
%------------------------------------------------------------------------------
function Flux = ComputeFluxes(alpha,beta,N,u,v,at_x,vel,scheme,dt,dx)
lambda=1.0;  

if scheme==0
   % Central-Diff, 
   temp = 0.5*(FLUX(u) + FLUX(v));
elseif scheme==1
   % Lax-friedrichs
   temp = 0.5*(FLUX(u) + FLUX(v)) - (0.5/lambda)*(v - u);
elseif scheme==2
   % Lax-wendroff
   a = ADFLUX(u,v);
   temp = 0.5*(FLUX(u) + FLUX(v)) - (0.5*lambda*a)*(FLUX(v) - FLUX(u));
elseif scheme==3
   % Murman-Roe
   a = ADFLUX(u,v);
   temp = 0.5*(FLUX(u) + FLUX(v)) - 0.5*abs(a)*(v - u);
elseif scheme==4
   % Godunov flux
   if u < v
      u1 = max(0,u);
      u2 = min(0,v);
      temp = max( FLUX(u1), FLUX(u2) );
   else
      temp = max( FLUX(u), FLUX(v) );
   end
elseif scheme==5
   % Roe with entropy fix
   delta = 0.5*max(0, abs(v-u));
   a = abs(ADFLUX(u,v));
   if a < delta
      a = delta;
   end
   temp = 0.5*(FLUX(u) + FLUX(v)) - 0.5*a*(v - u);
end
   
   Flux=zeros((N+1),1);
  for i=1:N+1  
   P = JacobiP(at_x, alpha, beta, i-1);   
      Flux(i,1) =temp*P ;      
  end   
end



