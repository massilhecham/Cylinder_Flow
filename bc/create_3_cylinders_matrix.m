function M = create_3_cylinders_matrix(n, D, config, gap_ratio)
% Crée une matrice nxn représentant un domaine contenant 1 ou 3 cylindres.
% - M(i,j) = 0 : extérieur
% - M(i,j) = 1, 2, 3 : bord des cylindres (numérotés de gauche à droite)
% - M(i,j) = 4 : intérieur d’un cylindre
%
% Entrées :
%   n         : taille du maillage (n x n)
%   D         : diamètre d’un cylindre (en unité normalisée)
%   config    : configuration des cylindres ('single', 'tandem', 'side by side', 'triangle', 'etage')
%   gap_ratio : espacement entre les centres des cylindres en multiple de D

% Rayon du cylindre
radius = D / 2;

% Espacement entre les centres
gap = gap_ratio * D;

% Grille de coordonnées normalisées (entre 0 et 1)
[X, Y] = meshgrid(linspace(0, 1, n), linspace(0, 1, n));

% Matrice résultat
M = zeros(n, n);

% Tolérance pour déterminer le bord
tol = 1 / n;

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
        centers = [0.5 - gap/2, 0.5;
                   0.5 + gap/2, 0.5 + h;
                   0.5 + gap/2, 0.5 - h];

    case 'etage'
        L = gap_ratio * D;
        d = L / sqrt(2);
        centers = [0.5 + d, 0.5 - d;
                   0.5,     0.5;
                   0.5 - d, 0.5 + d];

    otherwise
        error('Configuration non reconnue. Choisissez "single", "tandem", "side by side", "triangle" ou "etage".');
end

% Trier les centres de gauche à droite, puis de haut en bas
[~, idx] = sortrows(centers, [1, 2]);
centers = centers(idx, :);

% Boucle sur chaque cylindre
for k = 1:size(centers, 1)
    cx = centers(k, 1);
    cy = centers(k, 2);

    % Boucle sur chaque point du domaine
    for i = 1:n
        for j = 1:n
            % Distance entre le point (i,j) et le centre (cx, cy)
            dist = sqrt((X(i, j) - cx)^2 + (Y(i, j) - cy)^2);

            if dist <= radius + tol && dist > radius
                M(i, j) = k;  % Bord du cylindre
            elseif dist <= radius
                M(i, j) = 4;  % Intérieur du cylindre
            end
        end
    end
end
end
