clear; close all; clc;

addpath basic
addpath bc
addpath turbulence
addpath Calculs

% === PARAMÈTRES D'ENTRÉE ===
U_p = 0.01629;                 % Vitesse du fluide en entrée
nu_p = 1.9548e-05;             % Viscosité cinématique physique
rho0 = 1;                      % Densité initiale
Diameter = 0.12;               % Diamètre du cylindre
gap_ratio0 = 1.5;              % Rapport d'espacement entre les cylindres
config = 'tandem';             % Configuration: 'single', 'tandem', 'side by side', 'triangle', 'etage'
nodes = 834;                   % Nombre de nœuds du maillage (maillage carré nodes x nodes)
dt = 0.00516;                  % Pas de temps
timesteps = 500000;           % Nombre total d'itérations
err = 0.0;                     % Critère d'arrêt (variation relative de la vitesse)
nutilde0 = 1e-5;               % Valeur initiale pour le champ de turbulence
dh = 1/(nodes-1);              % Espacement entre les nœuds

% === PARAMÈTRES DÉRIVÉS ===
Re = Diameter * U_p / nu_p;
disp(['Reynolds number: ' num2str(Re)]);
u_lb = U_p * dt / dh;
disp(['Lattice speed: ' num2str(u_lb)]);
nu_lb = nu_p * dt / dh^2;
disp(['Lattice viscosity: ' num2str(nu_lb)]);
tau = 3 * nu_lb + 0.5;
disp(['Original relaxation time: ' num2str(tau)]);
omega = 1 / tau;
disp(['Physical relaxation parameter: ' num2str(omega)]);
Re_lb = u_lb * (Diameter/dh) / nu_lb;
disp(['Reynolds number (lattice units): ' num2str(Re_lb)]);
time_factor = U_p * dt * timesteps / Diameter;
disp(['time_factor: ' num2str(time_factor)]);

% === INITIALISATION ===
rho = rho0 * ones(nodes, nodes);
u = zeros(nodes, nodes); v = zeros(nodes, nodes);
gap_ratio = gap_ratio0 + 1;
u(end, 2:end-1) = u_lb;

% Distributions initiales (équilibre)
f = compute_feq(rho, u, v);

% Conditions aux limites initiales
f = moving_wall_bc(f, 'north');
f = moving_wall_bc(f, 'south');
f = inlet_bc(f, u_lb, 'west');
f = outlet_bc(f, 'east');
f = cylinder_bc(f, Diameter, nodes, config, gap_ratio);

% Turbulence : distance aux murs et champ initial
d = compute_moving_wall_distances(nodes);
nutilde = nutilde0 * ones(nodes, nodes);
[omega, nut, nutilde] = update_nut(nutilde, nu_lb, dt, dh, d, u, v);

% Préparation des résultats
uu_prev = zeros(nodes, nodes);
cylinder = create_3_cylinders_matrix(nodes, Diameter, config, gap_ratio);
time = [];
CD_1_vals = []; CL_1_vals = [];
CD_2_vals = []; CL_2_vals = [];
CD_3_vals = []; CL_3_vals = [];

% === BOUCLE PRINCIPALE ===
disp(['Running ' num2str(timesteps) ' timesteps...']);
tic;
for iter = 1:timesteps
    if mod(iter, timesteps/10) == 0
        disp(['Ran ' num2str(iter) ' iterations']);
    end

    % Étape de collision (avec turbulence)
    f = collide_sa(f, u, v, rho, omega);

    % Conditions aux limites (pré-streaming)
    f = moving_wall_bc(f, 'north');
    f = moving_wall_bc(f, 'south');
    f = outlet_bc(f, 'east');
    f = inlet_bc(f, u_lb, 'west');
    f = cylinder_bc(f, Diameter, nodes, config, gap_ratio);

    % Étape de propagation (streaming)
    f = stream(f);

    % Conditions aux limites (post-streaming)
    f = moving_wall_bc(f, 'north');
    f = moving_wall_bc(f, 'south');
    f = outlet_bc(f, 'east');
    f = inlet_bc(f, u_lb, 'west');
    f = cylinder_bc(f, Diameter, nodes, config, gap_ratio);

    % Reconstruction des variables macroscopiques
    [u, v, rho] = reconstruct_macro_all(f);

    % Imposition des conditions sur les parois et l'entrée/sortie
    u([1 end], :) = 0;
    v([1 end], :) = 0;
    u(2:end-1,1) = u_lb; v(2:end-1,1) = 0;
    u(:, end) = u(:, end-1); v(:, end) = v(:, end-1);
    u(cylinder ~= 0) = 0; v(cylinder ~= 0) = 0;

    % Mise à jour du modèle de turbulence
    [omega, nut, nutilde] = update_nut(nutilde, nu_lb, dt, dh, d, u, v);

    % Calcul des forces aérodynamiques
    Coeffs = compute_forces_coeffs(f, nodes, Diameter, rho0, u_lb, cylinder);
    time = [time, iter * dt * U_p / Diameter];
    CD_1_vals = [CD_1_vals, Coeffs(1, 1)];
    CL_1_vals = [CL_1_vals, Coeffs(1, 2)];
    CD_2_vals = [CD_2_vals, Coeffs(2, 1)];
    CL_2_vals = [CL_2_vals, Coeffs(2, 2)];
    CD_3_vals = [CD_3_vals, Coeffs(3, 1)];
    CL_3_vals = [CL_3_vals, Coeffs(3, 2)];

    % === VISUALISATION INTERACTIVE ===
    if mod(iter, 10) == 0
        % Graphiques des coefficients
        figure(1); clf;
        plot(time, CD_1_vals, 'r-', time, CD_2_vals, 'g-', time, CD_3_vals, 'b-');
        xlabel("tU / D"); ylabel("CD"); legend("CD_1", "CD_2", "CD_3");
        title(sprintf("Coefficients de traînée en fonction du temps \n avec Re = %.0f, configuration %s, gap ratio = %.1f", Re, config, gap_ratio0));
        drawnow;

        figure(2); clf;
        plot(time, CL_1_vals, 'r--', time, CL_2_vals, 'g-', time, CL_3_vals, 'b-');
        xlabel("tU / D"); ylabel("CL"); legend("CL_1", "CL_2", "CL_3");
        title(sprintf("Coefficients de portnance en fonction du temps \n avec Re = %.0f, configuration %s, gap ratio = %.1f", Re, config, gap_ratio0));
        drawnow;

        % Affichage du champ de vitesse normalisé
        uu = sqrt(u.^2 + v.^2) / u_lb;
        uu(cylinder ~= 0) = NaN;
        figure(3); clf; imagesc(uu); colorbar; axis equal off;
        cb = colorbar;        % Ajoute le colorbar
        % Titre et labels des axes
        title(sprintf('Module de vitesse avec Re = %.0f, configuration %s, gap ratio = %.1f', Re, config, gap_ratio0));
        ylabel(cb, 'u / u_{lb}'); 
        drawnow;
    end

    % Vérification de la convergence
    uu = sqrt(u.^2 + v.^2) / u_lb;
    diff_max = max(max(100 * abs(uu - uu_prev) ./ (uu_prev + 1e-8), [], 'omitnan'));
    if diff_max <= err
        break;
    end
    uu_prev = uu;
end

% Calcul du champ de vorticité
w = vorticity(u, v, nodes);

% Temps d'exécution
elapsedTime = toc;
disp("Done!");
fprintf("Temps d'exécution total : %.4f secondes\n", elapsedTime);
