function f = cylinder_bc(f, D,nodes,config,gap_ratio)
  % Boucle sur tous les nœuds de la matrice
    M = create_3_cylinders_matrix(nodes, D, config, gap_ratio);

    for i = 2:nodes-1
        for j = 2:nodes-1
            if M(i, j) == 1  || M(i, j) == 2 || M(i, j) == 3% Si le point appartient au cercle
                % Bounce-back : inversion des distributions
                 f(i, j, 2) = f(i, j, 4);
                 f(i, j, 3) = f(i, j, 5);
                 f(i, j, 4) = f(i, j, 2);
                 f(i, j, 5) = f(i, j, 3);
                 f(i, j, 6) = f(i, j, 8);
                 f(i, j, 7) = f(i, j, 9);
                 f(i, j, 8) = f(i, j, 6);
                 f(i, j, 9) = f(i, j, 7);  

            end
        end
    end
end


