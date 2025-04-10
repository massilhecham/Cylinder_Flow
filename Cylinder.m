% UNDER CONSTRUCTION

% A Lattice Boltzmann (single relaxation time) D2Q9 solver,
% with the Spalart Allmaras turbulence model
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

%L_p = 4; %1.1; % Cavity dimension


tic;
% Paramètres d'entrée.

U_p =0.01225; %1.1; % Cavity lid velocity.
nu_p = 9.1875e-6; % 1.586e-5; % Physical kinematic viscosity.
rho0 = 1;
Diameter=0.075; % Diamètre du cylindre
gap_ratio = 1.5; % Espacement relatif des cylindres
config  = 'tandem';   % Configurations:
% 'single' (1 cylindre);
% 'tandem' (3 cylindres alignés horizontalement);
% 'side by side' (3 cylindres alignés verticalement);
% 'triangle' (3 cylindres formant un triangle équilatéral face à l'écoulement);

nodes = 1143;  % maillage du domaine (nodes x nodes)
dt = 0.005; % Pas de temps
timesteps = 400000; % nombres d'itérations sur le pas de temps
err = 0.0; % condition d'arrêt sur la boucle (écart relatif sur l'amplitude de la vitesse)

nutilde0 = 1e-5; % initial nutilde value (should be non-zero for seeding).
dh = 1/(nodes-1);

% Derived nondimensional parameters.
Re = (Diameter) * U_p / nu_p;
disp(['Reynolds number: ' num2str(Re)]);
% Derived physical parameters.
t_p = Diameter / U_p;
disp(['Physical time scale: ' num2str(t_p) ' s']);
% Derived discrete parameters.
u_lb = U_p*dt / dh;
disp(['Lattice speed: ' num2str(u_lb)])
% nu_lb = dt / dh^2 / Re;
nu_lb = nu_p*dt/(dh^2);
disp(['Lattice viscosity: ' num2str(nu_lb)]);
tau = 3*nu_lb + 0.5;
disp(['Original relaxation time: ' num2str(tau)]);
omega = 1 / tau;
disp(['Physical relaxation parameter: ' num2str(omega)]);
Re_lb = (u_lb)*(Diameter/dh)/nu_lb;
disp(['Reynolds number (lattice units): ' num2str(Re_lb)]);

time_factor = U_p*dt*timesteps/Diameter;

disp(['time_factor: ' num2str(time_factor)]);

% % Derived nondimensional parameters.
% Re = Diameter * U_p / nu_p;
% disp(['Reynolds number: ' num2str(Re)]);
% % Derived physical parameters.
% t_p = Diameter / U_p;
% disp(['Physical time scale: ' num2str(t_p) ' s']);
% % Derived discrete parameters.
% dh = 1/(nodes-1);
% nu_lb = Diameter*dt / dh^2 / Re;
% disp(['Lattice viscosity: ' num2str(nu_lb)]);
% tau = 3*nu_lb + 0.5;
% disp(['Original relaxation time: ' num2str(tau)]);
% omega = 1 / tau;
% disp(['Physical relaxation parameter: ' num2str(omega)]);
% u_lb = dt / dh;
% disp(['Lattice speed: ' num2str(u_lb)])
% Re_lb = (u_lb)*(Diameter/dh)/nu_lb;
% disp(['Reynolds number (lattice units): ' num2str(Re_lb)]);

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
f = inlet_bc(f,u_lb,'west');
f = outlet_bc(f,'east');
f = cylinder_bc(f, Diameter,nodes,config,gap_ratio);
% Initialize turbulence stuff.
d = compute_wall_distances(nodes);
nutilde = nutilde0*ones(nodes,nodes);
[omega, nut, nutilde] = update_nut(nutilde,nu_lb,dt,dh,d,u,v);

uu_prev=zeros(nodes,nodes);


time = []; % Vecteur temps
CD_1_vals = [];
CL_1_vals = [];
CD_2_vals = [];
CL_2_vals = [];
CD_3_vals = [];
CL_3_vals = [];

cylinder = create_3_cylinders_matrix(nodes,Diameter,config,gap_ratio);

% Main loop.
disp(['Running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps
    if (mod(iter,timesteps/10)==0)
        disp(['Ran ' num2str(iter) ' iterations']);
    end

    % Collision.
    f = collide_sa(f, u, v, rho, omega); %collision

    % Apply meso BCs.
    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = outlet_bc(f,'east');
    f = inlet_bc(f,u_lb,'west');
    f = cylinder_bc(f, Diameter, nodes,config,gap_ratio);
    % Streaming.
    f = stream(f);

    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = outlet_bc(f,'east');
    f = inlet_bc(f,u_lb,'west');
    f = cylinder_bc(f, Diameter,nodes,config,gap_ratio);

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

    % Sortie à droite (Neumann)
    u(:,end) = u(:,end-1);  % Condition de sortie (extrapolation)
    v(:,end) = v(:,end-1);  % Condition de sortie (extrapolation)





    for i = 1:nodes
        for j = 1:nodes
            if cylinder(i,j) == 1 || cylinder(i,j) == 2 || cylinder(i,j) == 3 || cylinder(i,j) == 4
                u(i,j) = 0;
                v(i,j) = 0;
            end
        end
    end

    [omega, nut, nutilde] = update_nut(nutilde,nu_lb,dt,dh,d,u,v);


    Coeffs = compute_forces_coeffs(f, nodes, Diameter, rho0, u_lb, config, gap_ratio);
    time = [time, iter*dt*U_p/Diameter;];
    CD_1_vals = [CD_1_vals, Coeffs(1, 1)];
    CL_1_vals = [CL_1_vals, Coeffs(1, 2)];
    CD_2_vals = [CD_2_vals, Coeffs(2, 1)];
    CL_2_vals = [CL_2_vals, Coeffs(2, 2)];
    CD_3_vals = [CD_3_vals, Coeffs(3, 1)];
    CL_3_vals = [CL_3_vals, Coeffs(3, 2)];

    % VISUALIZATION
    % Modified from Jonas Latt's cavity code on the Palabos website.
    if (mod(iter, 10) == 0)

        

        % Figure pour CD
        figure(1);
        clf;
        hold on;
        plot(time, CD_1_vals, 'r-', 'DisplayName', 'CD_1');
        plot(time, CD_2_vals, 'g-', 'DisplayName', 'CD_2');
        plot(time, CD_3_vals, 'b-', 'DisplayName', 'CD_3');
        xlabel("tU / D");
        ylabel('Coefficients de traînée (CD)');
       % ylim([-5e-3 5e-3]);  % Limite de l'axe y
        legend;
        title(sprintf("Coefficients de traînée en fonction du temps \n avec Re = %.0f, configuration %s, gap ratio = %.1f", Re, config, gap_ratio));
        hold off;
        drawnow;

        % Figure pour CL
        figure(2);
        clf;
        hold on;
        plot(time, CL_1_vals, 'r-', 'DisplayName', 'CL_1');
        plot(time, CL_2_vals, 'g-', 'DisplayName', 'CL_2');
        plot(time, CL_3_vals, 'b-', 'DisplayName', 'CL_3');
        xlabel("tU / D");
        ylabel('Coefficients de portance (CL)');
       % ylim([-2e-3 2e-3]);  % Limite de l'axe y
        legend;
        title(sprintf("Coefficients de portance en fonction du temps \n avec Re = %.0f, configuration %s, gap ratio = %.1f", Re, config, gap_ratio));
        hold off;
        drawnow;


        uu = sqrt(u.^2 + v.^2) / u_lb;

        uu(cylinder ~= 0) = NaN;
        figure(3);
        imagesc(uu);

        colorbar;
        axis equal off;
        % colormap(jet);
        % set(gca, 'Color', [1 1 1]);
        drawnow;
        hold on;

        % Discrétisation du cercle
        theta = linspace(0, 2*pi, 100);

        % Calcul des valeurs nécessaires
        R = Diameter / dh / 2; % Rayon des cylindres en unités de grille
        gap = gap_ratio * Diameter / dh; % Espacement entre les cylindres en unités de grille
        xc = (nodes - 1) / 2; % Centre du domaine en x
        yc = (nodes - 1) / 2; % Centre du domaine en y

        switch config
            case 'single'
                x_cyl = xc;
                y_cyl = yc;

            case 'tandem'
                x_cyl = [xc - gap, xc, xc + gap];
                y_cyl = [yc, yc, yc];

            case 'side by side'
                x_cyl = [xc, xc, xc];
                y_cyl = [yc + gap, yc, yc - gap];

            case 'triangle'  % Triangle équilatéral pointant vers la gauche (inlet)
                h = sqrt(3) * (Diameter / dh / 2);  % Hauteur du triangle équilatéral
                % Centres des trois cylindres pour former un triangle équilatéral pointant vers la gauche
                x_cyl = [xc - gap/2 , xc + gap/2, xc + gap/2];  % Coordonnée x des 3 cylindres
                y_cyl = [yc , yc + h, yc - h];  % Coordonnée y des 3 cylindres

            otherwise
                error('Configuration non reconnue.');
        end

        % Colorier les cylindres en blanc
        for i = 1:length(x_cyl)
            % Dessiner chaque cylindre avec la bonne position
            fill(x_cyl(i) + R*cos(theta), y_cyl(i) + R*sin(theta), 'w', 'EdgeColor', 'none');
        end

        hold off;





    end



    uu = sqrt(u.^2+v.^2) / u_lb;
    uu_act= uu;
    epsilon = 1e-8;
    diff_max = max(max(100 * abs(uu_act - uu_prev) ./ (uu_prev + epsilon), [], 'omitnan'));

    if diff_max<= err
        break
    end
    uu_prev=uu_act;


end



% % Affichage du graphique des coefficients
% figure(2);
% clf;
% hold on;
% plot(time, CD_1_vals, 'r-', 'DisplayName', 'CD_1');
% plot(time, CL_1_vals, 'r--', 'DisplayName', 'CL_1');
% plot(time, CD_2_vals, 'g-', 'DisplayName', 'CD_2');
% plot(time, CL_2_vals, 'g--', 'DisplayName', 'CL_2');
% plot(time, CD_3_vals, 'b-', 'DisplayName', 'CD_3');
% plot(time, CL_3_vals, 'b--', 'DisplayName', 'CL_3');
% xlabel("tU / D");
% ylabel('Valeurs des coefficients');
% legend;
% title(sprintf("Coefficients aérodynamiques en fonction du nombre d'itérations \n avec Re = %.0f, une configuration %s et un rapport d'espacement de %.1f", Re, config,gap_ratio));
% drawnow;





pente = diff(CL_1_vals) ./ diff(time*dt);
regression = polyfit(time(500:end)*dt,CL_1_vals(500:end),1);
w = vorticity(u, v, nodes);



elapsedTime = toc;
disp('Done!');
fprintf('Temps d''exécution total : %.4f secondes\n', elapsedTime);

