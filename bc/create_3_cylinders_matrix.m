function M = create_3_cylinders_matrix(n, D, config, gap_ratio)
    % Crée une matrice nxn avec différentes configurations de cylindres.
    % Bord du cylindre = 1, 2 ou 3 selon le cylindre, intérieur du cylindre = 4, extérieur = 0.
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
    
    % Tolérance pour identifier le bord
    tol = 1 / n;

    % Définition des centres des cylindres en fonction de la configuration
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

        case 'triangle'  % Triangle équilatéral pointant vers la gauche (inlet)
            h = sqrt(3) * D / 2;  % Hauteur du triangle équilatéral
            % Centres des trois cylindres pour former un triangle équilatéral pointant vers la gauche
            centers = [0.5 - gap / 2, 0.5;                     % Premier cylindre à la pointe du triangle (à gauche)
                       0.5 + gap / 2, 0.5 + h ;             % Deuxième cylindre en haut à droite
                       0.5 + gap / 2, 0.5 - h ];            % Troisième cylindre en bas à droite

        otherwise
            error('Configuration non reconnue. Choisissez "single", "tandem", "side_by_side" ou "triangle".');
    end

    % Trier les centres des cylindres d'abord par leur coordonnée x (de gauche à droite)
    % En cas d'égalité, trier par la coordonnée y (de haut en bas)
    [~, idx] = sortrows(centers, [1, 2]);  % Trier d'abord par x, puis par y
    centers = centers(idx, :);

    % Numérotation des cylindres de gauche à droite, puis de haut en bas
    cyl_num = 1;  % Compteur pour la numérotation des cylindres
    for k = 1:size(centers, 1)
        cx = centers(k, 1);
        cy = centers(k, 2);

        % Parcours de la matrice pour définir les bords et l'intérieur des cylindres
        for i = 1:n
            for j = 1:n
                % Calculer la distance par rapport au centre du cylindre
                dist = sqrt((X(i, j) - cx)^2 + (Y(i, j) - cy)^2);

                % Assigner les valeurs aux cellules en fonction de la distance
                if dist <= radius + tol && dist > radius  % Bord du cylindre
                    M(i, j) = cyl_num;  % Bord du cylindre selon l'ordre de numérotation
                elseif dist <= radius  % Intérieur du cylindre
                    M(i, j) = 4;  % Intérieur du cylindre
                end
            end
        end

        % Incrémenter le compteur de numérotation
        cyl_num = cyl_num + 1;
    end
end



