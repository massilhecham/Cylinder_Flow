function Coeffs = compute_forces_coeffs(f, dt, nodes, D, rho, U_p, u_lb, config, gap_ratio)
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
                is_outward = check_outward_direction(cylinder, i, j, c(beta, :),D,gap_ratio,config);

                if is_outward == 1
                    % Vérifier que la norme n'est pas nulle pour éviter la division par zéro
                    if norm(c(beta, :)) ~= 0
                        opposite_beta = opposite_direction(beta);

                        % Calcul des distributions incidentes et rebondies
                        momentum_incident_D_1 = f(i, j, beta) * c(beta, 1);  % Distribution incidente
                        momentum_rebound_D_1 = f(i, j, opposite_beta) * c(opposite_beta, 1);  % Distribution après rebond

                        momentum_incident_L_1 = f(i, j, beta) * c(beta, 2);  % Distribution incidente
                        momentum_rebound_L_1 = f(i, j, opposite_beta) * c(opposite_beta, 2);  % Distribution après rebond

                        % Échange de quantité de mouvement
                        delta_momentum_D_1 = (momentum_rebound_D_1 - momentum_incident_D_1) / norm(c(beta, :));
                        delta_momentum_L_1 = (momentum_rebound_L_1 - momentum_incident_L_1) / norm(c(beta, :));

                        % Calcul des composantes de force pour le cylindre concerné
                        Force_Drag_1 = Force_Drag_1 + delta_momentum_D_1;
                        Force_Lift_1 = Force_Lift_1 + delta_momentum_L_1;

                    end
                end
            end
        end
       if cylinder(i, j) == 2
            for beta = 2:9
                is_outward = check_outward_direction(cylinder, i, j, c(beta, :),D,gap_ratio,config);

                if is_outward == 1
                    % Vérifier que la norme n'est pas nulle pour éviter la division par zéro
                    if norm(c(beta, :)) ~= 0
                        opposite_beta = opposite_direction(beta);

                        % Calcul des distributions incidentes et rebondies
                        momentum_incident_D_2 = f(i, j, beta) * c(beta, 1);  % Distribution incidente
                        momentum_rebound_D_2 = f(i, j, opposite_beta) * c(opposite_beta, 1);  % Distribution après rebond

                        momentum_incident_L_2 = f(i, j, beta) * c(beta, 2);  % Distribution incidente
                        momentum_rebound_L_2 = f(i, j, opposite_beta) * c(opposite_beta, 2);  % Distribution après rebond

                        % Échange de quantité de mouvement
                        delta_momentum_D_2 = (momentum_rebound_D_2 - momentum_incident_D_2) / norm(c(beta, :));
                        delta_momentum_L_2 = (momentum_rebound_L_2 - momentum_incident_L_2) / norm(c(beta, :));

                        % Calcul des composantes de force pour le cylindre concerné
                        Force_Drag_2 = Force_Drag_2 + delta_momentum_D_2;
                        Force_Lift_2 = Force_Lift_2 + delta_momentum_L_2;

                    end
                end
            end
       end
       if cylinder(i, j) == 3
            for beta = 2:9
                is_outward = check_outward_direction(cylinder, i, j, c(beta, :),D,gap_ratio,config);

                if is_outward == 1
                    % Vérifier que la norme n'est pas nulle pour éviter la division par zéro
                    if norm(c(beta, :)) ~= 0
                        opposite_beta = opposite_direction(beta);

                        % Calcul des distributions incidentes et rebondies
                        momentum_incident_D_3 = f(i, j, beta) * c(beta, 1);  % Distribution incidente
                        momentum_rebound_D_3 = f(i, j, opposite_beta) * c(opposite_beta, 1);  % Distribution après rebond

                        momentum_incident_L_3 = f(i, j, beta) * c(beta, 2);  % Distribution incidente
                        momentum_rebound_L_3 = f(i, j, opposite_beta) * c(opposite_beta, 2);  % Distribution après rebond

                        % Échange de quantité de mouvement
                        delta_momentum_D_3 = (momentum_rebound_D_3 - momentum_incident_D_3) / norm(c(beta, :));
                        delta_momentum_L_3 = (momentum_rebound_L_3 - momentum_incident_L_3) / norm(c(beta, :));

                        % Calcul des composantes de force pour le cylindre concerné
                        Force_Drag_3 = Force_Drag_3 + delta_momentum_D_3;
                        Force_Lift_3 = Force_Lift_3 + delta_momentum_L_3;

                    end
                end
            end
        end
    end
end

% Mise à l'échelle des forces
Force_Drag_1 = Force_Drag_1 * (dx^2 / dt);
Force_Lift_1 = Force_Lift_1 * (dx^2 / dt);
Force_Drag_2 = Force_Drag_2 * (dx^2 / dt);
Force_Lift_2 = Force_Lift_2 * (dx^2 / dt);
Force_Drag_3 = Force_Drag_3 * (dx^2 / dt);
Force_Lift_3 = Force_Lift_3 * (dx^2 / dt);

% Calcul des coefficients aérodynamiques pour chaque cylindre

CD_1 = Force_Drag_1/(0.5*rho*((u_lb)^2)*(D/dx));
CL_1 = Force_Lift_1/(0.5*rho*((u_lb)^2)*(D/dx));
CD_2 = Force_Drag_2/(0.5*rho*((u_lb)^2)*(D/dx)); %((U_p/u_lb)^2)*(D*dx)) avant! multiplier par (u_lb^4)/((U_p*dx)^2)
CL_2 = Force_Lift_2/(0.5*rho*((u_lb)^2)*(D/dx)); %((U_p/u_lb)^2)*(D*dx)) avant! multiplier par (u_lb^4)/((U_p*dx)^2)
CD_3 = Force_Drag_3/(0.5*rho*((u_lb)^2)*(D/dx)); %((U_p/u_lb)^2)*(D*dx)) avant! multiplier par (u_lb^4)/((U_p*dx)^2)
CL_3 = Force_Lift_3/(0.5*rho*((u_lb)^2)*(D/dx)); %((U_p/u_lb)^2)*(D*dx)) avant! multiplier par (u_lb^4)/((U_p*dx)^2)

Coeffs = [CD_1,CL_1;
          CD_2,CL_2;
          CD_3,CL_3];

end

