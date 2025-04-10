function Coeffs = compute_forces_coeffs(f, nodes, D, rho, u_lb, config, gap_ratio)
% Initialisation des forces
Force_Drag_1 = 0;
Force_Lift_1 = 0;
Force_Drag_2 = 0;
Force_Lift_2 = 0;
Force_Drag_3 = 0;
Force_Lift_3 = 0;


% Écart horizontal/vertical entre les noeuds
dx = 1/(nodes-1);



% Directions discrètes (vecteurs de vitesse) pour D2Q9
c = [
    0, 0;   % Stationnaire
    1, 0;   % Droite
    0, 1;   % Haut
    -1, 0;  % Gauche
    0, -1;  % Bas
    1, 1;   % Diagonale haut-droite
    -1, 1;  % Diagonale haut-gauche
    -1, -1; % Diagonale bas-gauche
    1, -1   % Diagonale bas-droite
    ];

% Correspondance des directions opposées en D2Q9
opposite_direction = [1, 4, 5, 2, 3, 8, 9, 6, 7];

% Créer la matrice des cylindres
cylinder = create_3_cylinders_matrix(nodes, D, config, gap_ratio);

% Parcourir les nœuds pour calculer les forces sur chaque cylindre
for i = 1:nodes
    for j = 1:nodes
        if cylinder(i, j) == 1
            for beta = 2:9
                opposite_beta = opposite_direction(beta);

                ni = i + c(beta,2);
                nj = j + c(beta,1);

                % Calcul des distributions incidentes et rebondies
                if cylinder(ni,nj) == 0

                    momentum_incident_1 = f(ni, nj, opposite_beta);
                    momentum_rebound_1 = f(i, j, beta);
                    % Échange de quantité de mouvement

                    delta_momentum_1 = momentum_rebound_1 - momentum_incident_1;
                    % Calcul des composantes de force pour le cylindre concerné
                    Force_Drag_1 = Force_Drag_1 + delta_momentum_1*c(beta,1);
                    Force_Lift_1 = Force_Lift_1 + delta_momentum_1*c(beta,2);
                end
            end
        end
        if cylinder(i, j) == 2
            for beta = 2:9

                opposite_beta = opposite_direction(beta);

                ni = i + c(beta,2);
                nj = j + c(beta,1);
                if cylinder(ni,nj) == 0
                    % Calcul des distributions incidentes et rebondies

                    momentum_incident_2 = f(ni,nj,opposite_beta);
                    momentum_rebound_2 = f(i,j,beta);

                    % Échange de quantité de mouvement
                    delta_momentum_2 = momentum_rebound_2 - momentum_incident_2;
                    % Calcul des composantes de force pour le cylindre concerné
                    Force_Drag_2 = Force_Drag_2 + delta_momentum_2 * c(beta,1);
                    Force_Lift_2 = Force_Lift_2 + delta_momentum_2 * c(beta,2);

                end
            end
        end
        if cylinder(i, j) == 3
            for beta = 2:9
                opposite_beta = opposite_direction(beta);

                ni = i + c(beta,2);
                nj = j + c(beta,1);
                if cylinder(ni,nj) == 0
                    % Calcul des distributions incidentes et rebondies

                    momentum_incident_3 = f(ni,nj,opposite_beta);
                    momentum_rebound_3 = f(i,j,beta);

                    % Échange de quantité de mouvement
                    delta_momentum_3 = momentum_rebound_3 - momentum_incident_3;
                    % Calcul des composantes de force pour le cylindre concerné
                    Force_Drag_3 = Force_Drag_3 + delta_momentum_3 * c(beta,1);
                    Force_Lift_3 = Force_Lift_3 + delta_momentum_3 * c(beta,2);
                end
            end
        end
    end
end

% Mise à l'échelle des forces

% Calcul des coefficients aérodynamiques pour chaque cylindre

CD_1 = Force_Drag_1/(0.5*rho*((u_lb)^2)*(D/dx));
CL_1 = Force_Lift_1/(0.5*rho*((u_lb)^2)*(D/dx));
CD_2 = Force_Drag_2/(0.5*rho*((u_lb)^2)*(D/dx));
CL_2 = Force_Lift_2/(0.5*rho*((u_lb)^2)*(D/dx));
CD_3 = Force_Drag_3/(0.5*rho*((u_lb)^2)*(D/dx));
CL_3 = Force_Lift_3/(0.5*rho*((u_lb)^2)*(D/dx));

Coeffs = [CD_1,CL_1;
          CD_2,CL_2;
          CD_3,CL_3];

end

