function is_outward = check_outward_direction(M, i, j, vitesse)
    % Dimensions de la matrice
    [n, ~] = size(M);
    
    % Définition du centre du domaine
    center_x = (n + 1) / 2;
    center_y = (n + 1) / 2;

    % Vérifier si le nœud est sur le bord du cercle
    if M(i, j) ~= 1
        is_outward = 0;
        return;
    end

    % Calcul de la normale extérieure (vers l'extérieur du cercle)
    normal = [i - center_x, j - center_y];

    % Produit scalaire entre la direction et la normale
    if dot(vitesse, normal) > 0
        is_outward = 1;
    else
        is_outward = 0;
    end
end
