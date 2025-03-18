function M = create_circle_matrix(n, D)
    % Crée une matrice nxn avec un cercle de diamètre D au centre
    % Bord du cercle = 1, intérieur du cercle = 2, extérieur = 0

    % Rayon du cercle
    radius = D / 2;  

    % Coordonnées normalisées [0,1]
    [X, Y] = meshgrid(linspace(0, 1, n), linspace(0, 1, n));

    % Centre du cercle (au milieu du domaine)
    center = 0.5;  

    % Distance radiale par rapport au centre
    r = sqrt((X - center).^2 + (Y - center).^2);

    % Tolérance pour identifier la bordure
    tol = 1.5 / n;  

    % Matrice de sortie
    M = zeros(n, n);

    % Marquer l'intérieur du cercle avec 2
    M(r < radius) = 2;

    % Marquer la bordure du cercle avec 1
    M(abs(r - radius) < tol) = 1;
end



