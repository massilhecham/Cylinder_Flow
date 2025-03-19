function [f,u,v] = cylinder_bc(f, D, u, v,nodes,config)
  % Boucle sur tous les nœuds de la matrice
    M = create_circle_matrix(nodes, D);

    for i = 2:nodes-1
        for j = 2:nodes-1
            if M(i, j) == 1  % Si le point appartient au cercle
                % Bounce-back : inversion des distributions
                f(i, j, 1) = f(i, j, 3);
                f(i, j, 2) = f(i, j, 4);
                f(i, j, 3) = f(i, j, 1);
                f(i, j, 4) = f(i, j, 2);
                f(i, j, 5) = f(i, j,7);
                f(i, j, 6) = f(i, j,8);
                f(i, j, 7) = f(i, j,5);
                f(i, j, 8) = f(i, j,6);
                f(i, j, 9) = f(i, j,9);  % (stationnaire, ne change pas)
            % elseif M(i, j) == 2 %Points situées dans le cercle
            %     u(i,j) = 0;
            %     v(i,j) = 0;
            end
        end
    end
end
