function w = vorticity(u, v, N)
    % Calcul de la vorticité w = dv/dx - du/dy
    % u : matrice des vitesses horizontales (N x N)
    % v : matrice des vitesses verticales (N x N)
    % N : nombre de points dans chaque direction (grille N x N)
    % Domaine carré de côté 1 avec maillage uniforme

    % Pas spatial (identique en x et y)
    dx = 1 / (N-1);
    dy = dx;

    % Initialisation de la matrice de vorticité
    w = zeros(N, N);

    % ====== Différences centrées à l'intérieur du domaine ======
    w(2:end-1, 2:end-1) = (v(2:end-1, 3:end) - v(2:end-1, 1:end-2)) / (2 * dx) ...
                        - (u(3:end, 2:end-1) - u(1:end-2, 2:end-1)) / (2 * dy);

    % ====== Traitement des bords ======
    % Bord gauche (x = 1) → Différence avant en x, arrière en y
    w(:, 1) = (v(:, 2) - v(:, 1)) / dx - (u([2:end 1], 1) - u(:, 1)) / dy; 

    % Bord droit (x = N) → Différence arrière en x, avant en y
    w(:, end) = (v(:, end) - v(:, end-1)) / dx - (u(:, end) - u([end 1:end-1], end)) / dy; 

    % Bord bas (y = 1) → Différence avant en y
    w(1, :) = (v(2, :) - v(1, :)) / dx - (u(1, :) - u(1, :)) / dy; % u_y = 0 ici

    % Bord haut (y = N) → Différence arrière en y
    w(end, :) = (v(end, :) - v(end-1, :)) / dx - (u(end, :) - u(end-1, :)) / dy;
end