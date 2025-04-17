function Coeffs = compute_forces_coeffs(f, nodes, D, rho, u_lb, cylinder)
% Calcule les coefficients de traînée (CD) et de portance (CL) pour 3 cylindres
% en utilisant la méthode d'échange de quantité de mouvement (momentum exchange).
% Inputs :
% - f : fonctions de distribution (3D : nodes x nodes x 9 directions)
% - nodes : nombre de points dans chaque direction
% - D : diamètre du cylindre (en unités physiques)
% - rho : densité du fluide
% - u_lb : vitesse de référence en unités lattice
% - cylinder : matrice d'identification des cylindres (1, 2, 3 pour les bords)

% Initialisation des forces pour les 3 cylindres
Force_Drag_1 = 0;  Force_Lift_1 = 0;
Force_Drag_2 = 0;  Force_Lift_2 = 0;
Force_Drag_3 = 0;  Force_Lift_3 = 0;

dx = 1 / (nodes - 1);  % Pas de grille

% Directions de D2Q9
c = [ 0,  0;
      1,  0;
      0,  1;
     -1,  0;
      0, -1;
      1,  1;
     -1,  1;
     -1, -1;
      1, -1 ];

opposite_direction = [1, 4, 5, 2, 3, 8, 9, 6, 7];

% Parcours de tous les nœuds du domaine
for i = 1:nodes
    for j = 1:nodes

        % Cylindre 1
        if cylinder(i, j) == 1
            for beta = 2:9
                opposite_beta = opposite_direction(beta);
                ni = i + c(beta,2);
                nj = j + c(beta,1);

                if cylinder(ni,nj) == 0  % Voisin en dehors du cylindre
                    momentum_incident_1 = f(ni, nj, opposite_beta);
                    momentum_rebound_1 = f(i, j, beta);
                    delta_momentum_1 = momentum_rebound_1 - momentum_incident_1;

                    Force_Drag_1 = Force_Drag_1 + delta_momentum_1 * c(beta,1);
                    Force_Lift_1 = Force_Lift_1 + delta_momentum_1 * c(beta,2);
                end
            end
        end

        % Cylindre 2
        if cylinder(i, j) == 2
            for beta = 2:9
                opposite_beta = opposite_direction(beta);
                ni = i + c(beta,2);
                nj = j + c(beta,1);

                if cylinder(ni,nj) == 0
                    momentum_incident_2 = f(ni,nj,opposite_beta);
                    momentum_rebound_2 = f(i,j,beta);
                    delta_momentum_2 = momentum_rebound_2 - momentum_incident_2;

                    Force_Drag_2 = Force_Drag_2 + delta_momentum_2 * c(beta,1);
                    Force_Lift_2 = Force_Lift_2 + delta_momentum_2 * c(beta,2);
                end
            end
        end

        % Cylindre 3
        if cylinder(i, j) == 3
            for beta = 2:9
                opposite_beta = opposite_direction(beta);
                ni = i + c(beta,2);
                nj = j + c(beta,1);

                if cylinder(ni,nj) == 0
                    momentum_incident_3 = f(ni,nj,opposite_beta);
                    momentum_rebound_3 = f(i,j,beta);
                    delta_momentum_3 = momentum_rebound_3 - momentum_incident_3;

                    Force_Drag_3 = Force_Drag_3 + delta_momentum_3 * c(beta,1);
                    Force_Lift_3 = Force_Lift_3 + delta_momentum_3 * c(beta,2);
                end
            end
        end
    end
end

% Mise à l'échelle des forces (formule standard des coefficients)
CD_1 = Force_Drag_1 / (0.5 * rho * u_lb^2 * (D / dx));
CL_1 = Force_Lift_1 / (0.5 * rho * u_lb^2 * (D / dx));
CD_2 = Force_Drag_2 / (0.5 * rho * u_lb^2 * (D / dx));
CL_2 = Force_Lift_2 / (0.5 * rho * u_lb^2 * (D / dx));
CD_3 = Force_Drag_3 / (0.5 * rho * u_lb^2 * (D / dx));
CL_3 = Force_Lift_3 / (0.5 * rho * u_lb^2 * (D / dx));

% Regroupement des coefficients dans une matrice 3x2
Coeffs = [CD_1, CL_1;
          CD_2, CL_2;
          CD_3, CL_3];
end
