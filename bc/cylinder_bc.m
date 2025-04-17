function f = cylinder_bc(f, D, nodes, config, gap_ratio)
% Applique la condition de rebond (bounce-back) sur les bords des cylindres
% dans le domaine de simulation. Le masque des cylindres est généré selon la
% configuration souhaitée (tandem, triangle, etc.).

% Génération de la matrice des cylindres (masque)
M = create_3_cylinders_matrix(nodes, D, config, gap_ratio);

% Directions discrètes (D2Q9)
c = [ 0,  0;   % 1 - au repos
      1,  0;   % 2 - droite
      0,  1;   % 3 - haut
     -1,  0;   % 4 - gauche
      0, -1;   % 5 - bas
      1,  1;   % 6 - haut droite
     -1,  1;   % 7 - haut gauche
     -1, -1;   % 8 - bas gauche
      1, -1 ]; % 9 - bas droite

% Direction opposée pour chaque direction (utile pour bounce-back)
opposite_direction = [1, 4, 5, 2, 3, 8, 9, 6, 7];

% Parcours des nœuds internes
for i = 2:nodes-1
    for j = 2:nodes-1
        % Vérifie si le nœud est un bord de cylindre
        if M(i, j) == 1 || M(i, j) == 2 || M(i, j) == 3
            for beta = 2:9
                % Vérifie si le vecteur vitesse sort du cylindre (pas bounce-back)
                is_outward = check_outward_direction(M, i, j, c(beta, :), D, gap_ratio, config);

                % Si la direction est entrante (vers l'intérieur du cylindre)
                if is_outward == 0
                    % Applique la condition bounce-back : réflexion dans la direction opposée
                    f(i, j, beta) = f(i, j, opposite_direction(beta));
                end
            end
        end
    end
end
end



