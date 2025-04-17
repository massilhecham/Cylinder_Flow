function is_outward = check_outward_direction(M, i, j, vitesse, D, gap_ratio, config)
    % Détermine si une direction donnée (vitesse) pointe vers l’extérieur d’un cylindre.
    % Sert à appliquer le rebond conditionnel dans les BC.
    % Entrées :
    % - M         : matrice des cylindres
    % - i, j      : position du point sur le bord du cylindre
    % - vitesse   : vecteur directionnel à tester (du schéma D2Q9)
    % - D         : diamètre du cylindre
    % - gap_ratio : rapport d’espacement entre les centres
    % - config    : disposition des cylindres ('single', 'tandem', etc.)

    % Taille de la matrice
    [n, ~] = size(M);
    
    % Calcul de l’espacement entre les centres
    gap = gap_ratio * D;
    
    % Définir les centres des cylindres selon la configuration choisie
    switch config
        case 'single'
            centers = [0.5, 0.5];

        case 'tandem'
            centers = [0.5 - gap, 0.5;
                       0.5,      0.5;
                       0.5 + gap, 0.5];

        case 'side by side'
            centers = [0.5, 0.5 - gap;
                       0.5, 0.5;
                       0.5, 0.5 + gap];

        case 'triangle'
            h = sqrt(3) * D / 2;
            centers = [0.5 - gap / 2, 0.5;
                       0.5 + gap / 2, 0.5 + h;
                       0.5 + gap / 2, 0.5 - h];

        case 'etage'
            L = gap_ratio * D;
            dir = L / sqrt(2);
            x2 = 0.5; y2 = 0.5;
            x1 = x2 + dir; y1 = y2 - dir;
            x3 = x2 - dir; y3 = y2 + dir;
            centers = [x1, y1;
                       x2, y2;
                       x3, y3];

        otherwise
            error('Configuration non reconnue. Choisissez une config valide.');
    end

    % Convertit les coordonnées physiques des centres en indices de matrice
    centers = round(centers * (n - 1) + 1);

    % Si le point (i,j) n’est pas un bord de cylindre, retour immédiat
    if M(i, j) ~= 1 && M(i, j) ~= 2 && M(i, j) ~= 3
        is_outward = 0;
        return;
    end

    % Identifier le centre de cylindre le plus proche
    dists = sqrt((centers(:,1) - j).^2 + (centers(:,2) - i).^2);
    [~, idx] = min(dists);
    center_x = centers(idx, 1);
    center_y = centers(idx, 2);

    % Calcul du vecteur normal pointant vers l’extérieur du cylindre
    normal = [j - center_x, i - center_y];
    normal = normal / norm(normal);

    % Produit scalaire pour vérifier si la vitesse pointe dans la même direction
    if norm(vitesse) > 0 && dot(vitesse, normal) > 0
        is_outward = 1; % Vitesse vers l’extérieur
    else
        is_outward = 0; % Vitesse vers l’intérieur ou tangentielle
    end
end
