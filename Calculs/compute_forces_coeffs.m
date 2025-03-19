function [Force_Lift, Force_Drag, Coeff_L, Coeff_D] = compute_forces_coeffs(f, dt, nodes, D, rho, U)
    % Calcule les forces et coefficients aérodynamiques sur un cylindre centré

    % Rayon et taille de grille
    R = D / 2;
    dx = 1 / (nodes - 1); % Espacement du maillage


    scaling_factor = dx^2 / dt;
    if scaling_factor > 1e6  % Seuil arbitraire à ajuster
        warning('Le facteur de mise à l’échelle est très grand, vérifiez dx et dt.');
    end
    
    c = zeros(9,2);
    c(1,:) = [0, 0];
    c(2,:) = [1, 0];
    c(3,:) = [0, 1];
    c(4,:) = [-1, 0];
    c(5,:) = [0, -1];
    c(6,:) = [1, 1];
    c(7,:) = [-1, 1];
    c(8,:) = [-1, -1];
    c(9,:) = [1, -1];


    % Création des coordonnées normalisées du domaine
    [X, Y] = meshgrid(linspace(0, 1, nodes), linspace(0, 1, nodes));

    % Centre du cylindre
    xc = 0.5;
    yc = 0.5;

    % Distance de chaque point au centre
    r = sqrt((X - xc).^2 + (Y - yc).^2);

    % **Correction ici**: on s'assure que `tol` est bien un SCALAIRE
    tol = dx * 1.5;

    % Détection des points de la surface du cylindre
    surface_mask = abs(r - R) < tol;

    % Sélection des indices correspondant à la surface du cylindre
    [I, J] = find(surface_mask);

    M=create_circle_matrix(nodes,D);

    Force_Lift=0;
    Force_Drag=0;

    for a=1:nodes
        for b=1:nodes
            
            if M(a,b)==1
                for beta = 1:9

                    if c(beta,:) == [0, 0]

                        Force_Drag= Force_Drag+ 2 * c(beta,1) * f(a, b, beta);

                    else
                        
                        Force_Drag= Force_Drag+ (2 * c(beta,1) * f(a, b, beta))/norm(c(beta,:));
                    end

                    if c(beta,2) == [0, 0]

                        Force_Lift= Force_Lift+ (2 * c(beta,2) * f(a, b, beta));

                    else

                        Force_Lift= Force_Lift+ (2 * c(beta,2) * f(a, b, beta))/norm(c(beta,:));
                    end
                end
            end
        end
    end
    
    Force_Drag = Force_Drag * (dx^2 / dt);
    Force_Lift = Force_Lift * (dx^2 / dt);

    Coeff_D = Force_Drag/(0.5*rho*(U^2)*D);
    Coeff_L = Force_Lift/(0.5*rho*(U^2)*D);

    % Initialisation des forces **comme scalaires**
    F_L = 0;
    F_D = 0;

    % Vérification de la taille de f
    [Nx, Ny, Q] = size(f);
    if Nx ~= nodes || Ny ~= nodes || Q ~= 9
        error('La taille de f est incorrecte ! Vérifiez que f est bien (nodes, nodes, 9).');
    end

    % Boucle sur les points de la surface du cylindre
    for k = 1:length(I)
        i = I(k);
        j = J(k);

        % **Forces locales sous forme de scalaires**
        force_local_L = 0; % Portance (direction y)
        force_local_D = 0; % Traînée (direction x)

        % Contribution de chaque direction de vitesse
        for alpha = 1:9
            force_local_D = force_local_D + 2 * c(alpha,1) * f(i, j, alpha);
            force_local_L = force_local_L + 2 * c(alpha,2) * f(i, j, alpha);
        end

        % **Ajout des forces sous forme SCALAIRE**
        F_L = F_L + (dx^2 / dt) * force_local_L;
        F_D = F_D + (dx^2 / dt) * force_local_D;
    end

    % Calcul des coefficients de portance et traînée
    dynamic_pressure = 0.5 * rho * U^2;
    
    if dynamic_pressure * D ~= 0
        C_L = F_L / (dynamic_pressure * D);  % Division SCALAIRE
        C_D = F_D / (dynamic_pressure * D);
    else
        C_L = NaN;
        C_D = NaN;
        warning('Attention : division par zéro lors du calcul de C_L et C_D');
    end
end
