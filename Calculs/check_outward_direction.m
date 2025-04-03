function is_outward = check_outward_direction(M, i, j, vitesse, D, gap_ratio, config)
    % Dimensions de la matrice
    [n, ~] = size(M);
    
    % Identifier les centres des cylindres en fonction de la configuration
    gap = gap_ratio * D;
    
    switch config
        case 'single'
            centers = [0.5, 0.5];
        
        case 'tandem'  % 3 cylindres alignés horizontalement
            centers = [0.5 - gap, 0.5;
                       0.5, 0.5;
                       0.5 + gap, 0.5];
        
        case 'side by side'  % 3 cylindres alignés verticalement
            centers = [0.5, 0.5 - gap;
                       0.5, 0.5;
                       0.5, 0.5 + gap];
        
        case 'triangle'  % Triangle équilatéral pointant vers la gauche
            h = sqrt(3) * D / 2;  % Hauteur du triangle
            centers = [0.5 - gap / 2, 0.5;
                       0.5 + (gap / 2), 0.5 + h;
                       0.5 + (gap / 2), 0.5 - h];
        
        otherwise
            error('Configuration non reconnue. Choisissez "single", "tandem", "side by side" ou "triangle".');
    end

    % Convertir les centres en indices de la grille
    centers = round(centers * (n - 1) + 1);  % Conversion en indices de la matrice (1-based)

    % Si la cellule (i, j) n'est pas sur le bord d'un cylindre (M(i,j) == 1, 2 ou 3)
    if M(i, j) ~= 1 && M(i, j) ~= 2 && M(i, j) ~= 3
        is_outward = 0;
        return;
    end
    
    % Trouver le centre du cylindre le plus proche
    dists = sqrt((centers(:,1) - j).^2 + (centers(:,2) - i).^2);  % Distance entre le point (i,j) et les centres
    [~, idx] = min(dists);  % Trouver l'indice du cylindre le plus proche
    center_x = centers(idx, 1);
    center_y = centers(idx, 2);
    
    % Calcul de la normale extérieure
    normal = [i - center_y, j - center_x];
    
    % Normalisation du vecteur normal
    normal = normal / norm(normal);
    
    % Vérifier si la vitesse pointe vers l'extérieur
    if norm(vitesse) > 0 && dot(vitesse, normal) > 0
        is_outward = 1;  % La vitesse pointe vers l'extérieur du cylindre
    else
        is_outward = 0;  % La vitesse ne pointe pas vers l'extérieur
    end
end



