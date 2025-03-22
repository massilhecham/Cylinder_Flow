
function is_outward = check_outward_direction(M, i, j, vitesse)
    % Dimensions de la matrice
    [n, ~] = size(M);
    
    % Définition du centre du domaine
    center_x = ceil(n / 2);
    center_y = ceil(n / 2);

    % Vérifier si le nœud est sur le bord du cylindre if M(i, j) ~= 1
    %     is_outward = 0; return;
    % end

    % Calcul de la normale extérieure (vectorisée)
    normal = [i - center_x, j - center_y];

    % Normalisation du vecteur normal
    normal = normal / norm(normal);

    % Vérifier si le vecteur vitesse pointe vers l'extérieur
    if norm(vitesse) > 0 && dot(vitesse, normal) >0
        is_outward = 1;
    else
        is_outward = 0;
    end
end
