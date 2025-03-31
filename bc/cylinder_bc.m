function f = cylinder_bc(f, D,nodes,config,gap_ratio)
% Boucle sur tous les nœuds de la matrice
M = create_3_cylinders_matrix(nodes, D, config, gap_ratio);


% Directions discrètes (vecteurs de vitesse) pour D2Q9
c = [
    0, 0;   % Stationnaire
    1, 0;   % Droite
    0, 1;   % Haut
    -1, 0;  % Gauche
    0, -1;  % Bas
    1, 1;   % Diagonale haut-droite
    -1, 1;  % Diagonale haut-gauche
    -1, -1; % Diagonale bas-gauche
    1, -1   % Diagonale bas-droite
    ];

opposite_direction = [1, 4, 5, 2, 3, 8, 9, 6, 7];

for i = 2:nodes-1
    for j = 2:nodes-1
        if M(i, j) == 1  || M(i, j) == 2 || M(i, j) == 3 % Si le point appartient au cercle

            for beta = 2:9

                is_outward = check_outward_direction(M, i, j, c(beta, :),D,gap_ratio,config);

                if is_outward == 0

                    f(i,j,beta) = f(i,j,opposite_direction(beta));

                    % Bounce-back : inversion des distributions
                    % f(i, j, 2) = f(i, j, 4);
                    % f(i, j, 3) = f(i, j, 5);
                    % f(i, j, 4) = f(i, j, 2);
                    % f(i, j, 5) = f(i, j, 3);
                    % f(i, j, 6) = f(i, j, 8);
                    % f(i, j, 7) = f(i, j, 9);
                    % f(i, j, 8) = f(i, j, 6);
                    % f(i, j, 9) = f(i, j, 7);
                end
            end
        end
    end
end
end


