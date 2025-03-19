% UNDER CONSTRUCTION

% A Lattice Boltzmann (single relaxation time) D2Q9 solver,
% with the Spalart Allmaras turbulence model, on a lid-driven cavity. 
% Cell centers (nodes) are placed on the boundaries. 
% Author: Robert Lee
% Email: rlee32@gatech.edu

clear;close all;clc;

addpath basic
addpath bc
addpath turbulence
addpath Calculs
% Algorithm steps:
% Initialize meso (f)
% Apply meso BCs
% Determine macro variables and apply macro BCs
% Loop:
%   Collide
%   Apply meso BCs
%   Stream
%   Apply meso BCs?
%   Determine macro variables and apply macro BCs
tic;
% Physical parameters.
L_p = 4; %1.1; % Cavity dimension. 
U_p = 1.2; %1.1; % Cavity lid velocity.
nu_p = 1.2e-3; % 1.586e-5; % Physical kinematic viscosity.
rho0 = 1;

Diameter=0.1; % Diamètre du cylindre


% Discrete/numerical parameters.
nodes =300;
dt = 0.0005;
timesteps = 300000;
nutilde0 = 1e-5; % initial nutilde value (should be non-zero for seeding).

% Derived nondimensional parameters.
Re = Diameter * U_p / nu_p;
disp(['Reynolds number: ' num2str(Re)]);
% Derived physical parameters.
t_p = L_p / U_p;
disp(['Physical time scale: ' num2str(t_p) ' s']);
% Derived discrete parameters.
dh = 1/(nodes-1);
nu_lb = dt / dh^2 / Re;
disp(['Lattice viscosity: ' num2str(nu_lb)]);
tau = 3*nu_lb + 0.5;
disp(['Original relaxation time: ' num2str(tau)]);
omega = 1 / tau;
disp(['Physical relaxation parameter: ' num2str(omega)]);
u_lb = dt / dh;
disp(['Lattice speed: ' num2str(u_lb)])

% if dt> (dh^2/(10*nu_lb))
%     error('dt> (dh^2/(10*nu_lb))');
% end
if nodes < Re/10
    error('nodes < Re/10');
end

% if (tau<0.5 | tau>2)
%     error('tau<0.5 | tau>2');
% end
% Determine macro variables and apply macro BCs
% Initialize macro, then meso.
rho = rho0*ones(nodes,nodes);
u = zeros(nodes,nodes);
v = zeros(nodes,nodes);

u(end,2:end-1) = u_lb;
% Initialize.
f = compute_feq(rho,u,v);
% Apply meso BCs.
f = wall_bc(f,'north');
f = wall_bc(f,'south');
f = outlet_bc(f,'east');
f = inlet_bc(f,u_lb,'west');
[f,u,v] = cylinder_bc(f, Diameter, u, v,nodes);
% Initialize turbulence stuff.
d = compute_wall_distances(nodes);
nutilde = nutilde0*ones(nodes,nodes);
[omega, nut, nutilde] = update_nut(nutilde,nu_lb,dt,dh,d,u,v);

uu_prev=zeros(nodes,nodes);
%[F_L, F_D, C_L, C_D] = compute_forces_coeffs(f, dt, nodes,Diameter, rho0, u_lb);
% Main loop.
disp(['Running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps
    if (mod(iter,timesteps/10)==0)
        disp(['Ran ' num2str(iter) ' iterations']);
    end
    
    % Collision.
    f = collide_sa(f, u, v, rho, omega);
    
    % Apply meso BCs.
    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = outlet_bc(f,'east');
    f = inlet_bc(f,u_lb,'west');
    [f,u,v] = cylinder_bc(f, Diameter, u, v,nodes);
    % Streaming.
    f = stream(f);
    
    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = outlet_bc(f,'east');
    f = inlet_bc(f,u_lb,'west');
    [f,u,v] = cylinder_bc(f, Diameter, u, v,nodes);
    
    % Determine macro variables and apply macro BCs
    [u,v,rho] = reconstruct_macro_all(f);
    % Murs en haut et en bas
    u(1,:) = 0;         % Mur en haut (u = 0)
    u(end,:) = 0;       % Mur en bas (u = 0)
    v(1,:) = 0;         % Mur en haut (v = 0)
    v(end,:) = 0;       % Mur en bas (v = 0)

    % Entrée à gauche
    u(2:end-1,1) = u_lb; % u_lb à l'entrée sauf aux coins
    v(2:end-1,1) = 0;    % Pas de composante verticale à l'entrée
    f = outlet_bc(f,'east');
    % Sortie à droite (Neumann)
    u(:,end) = u(:,end-1);  % Condition de sortie (extrapolation)
    v(:,end) = v(:,end-1);  % Condition de sortie (extrapolation)
    [f,u,v] = cylinder_bc(f, Diameter, u, v,nodes);
    cylinder = create_circle_matrix(nodes,Diameter);
    for i = 1:nodes
        for j = 1:nodes
            if cylinder(i,j) ==1
                u(i,j) = 0;
                v(i,j) = 0;
            end
        end
    end

    [omega, nut, nutilde] = update_nut(nutilde,nu_lb,dt,dh,d,u,v);

    % VISUALIZATION
    % Modified from Jonas Latt's cavity code on the Palabos website.
    if (mod(iter,10)==0)
        uu = sqrt(u.^2+v.^2) / u_lb;
         imagesc(flipud(uu));
%        imagesc(flipud(nut));
%         imagesc(flipud(omega));
        colorbar
        axis equal off; drawnow
    end
   uu = sqrt(u.^2+v.^2) / u_lb;
   uu_act= uu;
   epsilon = 1e-8;
   diff_max = max(max(100 * abs(uu_act - uu_prev) ./ (uu_prev + epsilon), [], 'omitnan'));

   if diff_max<=1
       break
   end
   uu_prev=uu_act;
end

w = vorticity(u, v, nodes);
[F_L, F_D, C_L, C_D] = compute_forces_coeffs(f, dt, nodes,Diameter, rho0, u_lb);
elapsedTime = toc;
disp('Done!');
fprintf('Temps d''exécution total : %.4f secondes\n', elapsedTime);
