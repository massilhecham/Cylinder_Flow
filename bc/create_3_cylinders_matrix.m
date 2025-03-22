function M = create_3_cylinders_matrix(n, D, config, gap_ratio)
    % Crée une matrice nxn avec différentes configurations de cylindres.
    % Bord du cylindre = 1, intérieur du cylindre = 2, extérieur = 0.
    %
    % Paramètres :
    % n : taille de la matrice (nxn)
    % D : diamètre des cylindres
    % config : type de configuration ('single', 'tandem', 'side_by_side', 'triangle')
    % gap_ratio : facteur de séparation entre les cylindres (ex: 1.5 signifie un espace de 0.5D)

    % Rayon du cylindre
    radius = D / 2;
    
    % Espacement entre les cylindres
    gap = gap_ratio * D;
    
    % Coordonnées normalisées [0,1]
    [X, Y] = meshgrid(linspace(0, 1, n), linspace(0, 1, n));

    % Matrice de sortie
    M = zeros(n, n);
    
    % Tolérance pour identifier la bordure
    tol = .5 / n;

    % Définition des centres des cylindres en fonction de la configuration
    switch config
        case 'single'
            centers = [0.5, 0.5];

        case 'tandem'  % 3 cylindres alignés horizontalement
            centers = [0.5 - gap, 0.5;
                       0.5, 0.5;
                       0.5 + gap, 0.5];

        case 'side_by_side'  % 3 cylindres alignés verticalement
            centers = [0.5, 0.5 - gap;
                       0.5, 0.5;
                       0.5, 0.5 + gap];

        case 'triangle'  % Triangle équilatéral pointant vers la gauche
            h = sqrt(3) * D / 2;  % Hauteur du triangle
            centers = [0.5 - gap, 0.5;
                       0.5 + (gap / 2), 0.5 - (h / 2);
                       0.5 + (gap / 2), 0.5 + (h / 2)];

        otherwise
            error('Configuration non reconnue. Choisissez "single", "tandem", "side_by_side" ou "triangle".');
    end

    % Boucle sur les centres des cylindres
    for k = 1:size(centers, 1)
        cx = centers(k, 1);
        cy = centers(k, 2);

        % Distance radiale par rapport au centre du cylindre
        r = sqrt((X - cx).^2 + (Y - cy).^2);

        % Marquer l'intérieur du cylindre avec 2
        M(r < radius) = 2;

        % Marquer la bordure du cylindre avec 1
        M(abs(r - radius) < tol) = 1;
    end
end
