function [Force_Lift, Force_Drag, Coeff_L, Coeff_D] = compute_forces_coeffs(f, dt, nodes, D, rho, U)
    % Calcule les forces et coefficients aérodynamiques sur un cylindre centré
    % Implémente la méthode d'échange de quantité de mouvement (Momentum Exchange Method)

    % Taille de grille
    dx = 1 / (nodes - 1); % Espacement du maillage

    % Vérification du facteur d'échelle
    scaling_factor = dx^2 / dt;
    if scaling_factor > 1e6  % Seuil arbitraire à ajuster
        warning('Le facteur de mise à l’échelle est très grand, vérifiez dx et dt.');
    end

    % Directions discrètes (vecteurs de vitesse) pour D2Q9
    c = [
        0,  0;  % Stationnaire
        1,  0;  % Droite
        0,  1;  % Haut
       -1,  0;  % Gauche
        0, -1;  % Bas
        1,  1;  % Diagonale haut-droite
       -1,  1;  % Diagonale haut-gauche
       -1, -1;  % Diagonale bas-gauche
        1, -1   % Diagonale bas-droite
    ];

    % Générer la matrice du cylindre
    M = create_circle_matrix(nodes, D);

    % Initialisation des forces
    Force_Lift = 0;
    Force_Drag = 0;

    % Parcourir les noeuds pour calculer les forces
    for a = 2:nodes-1
        for b = 2:nodes-1
            if M(a, b) == 1  % Bordure du cylindre

                for beta = 2:9
                    is_outward = check_outward_direction(M, a, b, c(beta,:));

                    if is_outward == 0

                        % Correspondance des directions opposées en D2Q9
                        opposite_direction = [1, 4, 5, 2, 3, 8, 9, 6, 7];
                        opposite_beta = opposite_direction(beta);

                        % Calcul des distributions incidentes et rebondies
                        momentum_incident_D = f(a, b, beta)*c(beta, 1);      % Distribution incidente
                        momentum_rebound_D = f(a, b, opposite_beta)*c(opposite_beta, 1);  % Distribution après rebond

                        momentum_incident_L = f(a, b, beta)*c(beta, 2);      % Distribution incidente
                        momentum_rebound_L = f(a, b, opposite_beta)*c(opposite_beta, 2);  % Distribution après rebond

                        % Échange de quantité de mouvement
                        delta_momentum_D = (momentum_rebound_D - momentum_incident_D)/norm(c(beta,:));
                        delta_momentum_L = (momentum_rebound_L - momentum_incident_L)/norm(c(beta,:));
                        % Calcul des composantes de force
                        Force_Drag = Force_Drag + delta_momentum_D ;
                        Force_Lift = Force_Lift + delta_momentum_L ;

                    end
                end
            end
        end
    end

    % Mise à l'échelle des forces
    Force_Drag = Force_Drag * (dx^2 / dt);
    Force_Lift = Force_Lift * (dx^2 / dt);

    % Calcul des coefficients aérodynamiques
    Coeff_D = Force_Drag / (0.5 * rho * (U^2) * D);
    Coeff_L = Force_Lift / (0.5 * rho * (U^2) * D);
end
