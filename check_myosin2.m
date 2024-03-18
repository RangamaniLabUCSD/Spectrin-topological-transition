function [r_s,myosin2,edges_s3,P] = check_myosin2(T_s,r_s,edges_s,edge_type,myosin2,edges_s1,edges_s3,P)
% check if myosin rods can attach to other spectrin triangles

% check the corresponding triangles
    aux_T1 = [find(T_s(P.myosin_T2,1) == edges_s(edges_s1,1));
        find(T_s(P.myosin_T2,2) == edges_s(edges_s1,1));
        find(T_s(P.myosin_T2,3) == edges_s(edges_s1,1))];
    aux_T2 = [find(T_s(P.myosin_T2,1) == edges_s(edges_s1,2));
        find(T_s(P.myosin_T2,2) == edges_s(edges_s1,2));
        find(T_s(P.myosin_T2,3) == edges_s(edges_s1,2))];
    aux_T3 = intersect(aux_T1,aux_T2);
    if ~isempty(aux_T3)
%             if the edge corresponded to the spectrin triangle, then check whether
%             the other end of the myosin rod can attach to another spectrin triangle
%             within the allowed lengths of myosin 
        aux_TT = [];
        aux_TT2 = [];
        for lll = 1:length(aux_T3)
            center_Tfree = (r_s(T_s(P.myosin_Tfree,1),:) + r_s(T_s(P.myosin_Tfree,2),:) +...
                r_s(T_s(P.myosin_Tfree,3),:))./3;
            aux = find(edge_type == 4);
            aux2 = [find(edges_s(aux,1) == myosin2(aux_T3(lll)));
                    find(edges_s(aux,2) == myosin2(aux_T3(lll)))];
            aux3 = [find(myosin2 == edges_s(aux(aux2),1));
                find(myosin2 == edges_s(aux(aux2),2))];
            aux4 = setdiff(aux3,aux_T3(lll));

            aux_d = center_Tfree - r_s(myosin2(aux4),:);
            aux_d = sqrt(dot(aux_d,aux_d,2));
            aux_new = find(aux_d < P.max_r2 & aux_d > P.min_r2);

            auxb = find(edge_type == 3);

            if ~isempty(aux_new)
%                     randomly select one of the spectrin traingles, trace the
%                     myosin rod and check if it intersects with the
%                     other myosins 
                aux_c = 1:length(aux_new);
                aux_r = randi(length(aux_c));

                a = (r_s(edges_s([aux;auxb],2),2) - r_s(edges_s([aux;auxb],1),2))./...
                    (r_s(edges_s([aux;auxb],2),1) - r_s(edges_s([aux;auxb],1),1));
                b = -a.*r_s(edges_s([aux;auxb],1),1) + r_s(edges_s([aux;auxb],1),2);
                r_s_aux = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
                c = (r_s(myosin2(aux4),2) - r_s_aux(2))./...
                    (r_s(myosin2(aux4),1) - r_s_aux(1));
                d = -c*r_s_aux(1) + r_s_aux(2);
                X = (b-d)./(c-a);
                Y = (b.*c - d.*a)./(c-a);
                aux_a = [find(a==inf);find(a==-inf)];
                if ~isempty(aux_a)
                    aux_aa = r_s(edges_s([aux;auxb],1),1);
                    X(aux_a) = aux_aa(aux_a);
                    Y(aux_a) = c*X(aux_a) + d;
                end
                if c == inf || c == -inf
                    X = r_s(myosin2(aux4),1)*ones(size(b));
                    Y = a.*X + b;
                end

                cond1 = (r_s(edges_s([aux;auxb],1),1)-P.d_inter < X).*(r_s(edges_s([aux;auxb],2),1)+P.d_inter > X);
                cond2 = (r_s(edges_s([aux;auxb],2),1)-P.d_inter < X).*(r_s(edges_s([aux;auxb],1),1)+P.d_inter > X);
                cond3 = (r_s_aux(1)-P.d_inter < X).*(r_s(myosin2(aux4),1) > X);
                cond4 = (r_s_aux(1)+P.d_inter > X).*(r_s(myosin2(aux4),1) < X);

                cond5 = (r_s(edges_s([aux;auxb],1),2)-P.d_inter < Y).*(r_s(edges_s([aux;auxb],2),2)+P.d_inter > Y);
                cond6 = (r_s(edges_s([aux;auxb],2),2)-P.d_inter < Y).*(r_s(edges_s([aux;auxb],1),2)+P.d_inter > Y);
                cond7 = (r_s_aux(2)-P.d_inter < Y).*(r_s(myosin2(aux4),2) > Y);
                cond8 = (r_s_aux(2)+P.d_inter > Y).*(r_s(myosin2(aux4),2) < Y);
                while (sum((cond1+cond2).*(cond3+cond4)) + sum((cond5+cond6).*(cond7+cond8)) > 0)...
                        && length(aux_c)>1 
                    
                    aux_c(aux_r) = [];
                    aux_r = randi(length(aux_c));

                    r_s_aux = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
                    c = (r_s(myosin2(aux4),2) - r_s_aux(2))./...
                        (r_s(myosin2(aux4),1) - r_s_aux(1));
                    d = -c*r_s_aux(1) + r_s_aux(2);

                    X = (b-d)./(c-a);
                    Y = (b.*c - d.*a)./(c-a);
                    aux_a = [find(a==inf);find(a==-inf)];
                    if ~isempty(aux_a)
                        aux_aa = r_s(edges_s([aux;auxb],1),1);
                        X(aux_a) = aux_aa(aux_a);
                        Y(aux_a) = c*X(aux_a) + d;
                    end
                    if c == inf || c == -inf
                        X = r_s(myosin2(aux4),1)*ones(size(b));
                        Y = a.*X + b;
                    end

                    cond1 = (r_s(edges_s([aux;auxb],1),1)-P.d_inter < X).*(r_s(edges_s([aux;auxb],2),1)+P.d_inter > X);
                    cond2 = (r_s(edges_s([aux;auxb],2),1)-P.d_inter < X).*(r_s(edges_s([aux;auxb],1),1)+P.d_inter > X);
                    cond3 = (r_s_aux(1)-P.d_inter < X).*(r_s(myosin2(aux4),1) > X);
                    cond4 = (r_s_aux(1)+P.d_inter > X).*(r_s(myosin2(aux4),1) < X);

                    cond5 = (r_s(edges_s([aux;auxb],1),2)-P.d_inter < Y).*(r_s(edges_s([aux;auxb],2),2)+P.d_inter > Y);
                    cond6 = (r_s(edges_s([aux;auxb],2),2)-P.d_inter < Y).*(r_s(edges_s([aux;auxb],1),2)+P.d_inter > Y);
                    cond7 = (r_s_aux(2)-P.d_inter < Y).*(r_s(myosin2(aux4),2) > Y);
                    cond8 = (r_s_aux(2)+P.d_inter > Y).*(r_s(myosin2(aux4),2) < Y);
                end
%                     check if the myosin rod does not intersect the other myosins and update its position, otherwise 
%                 remove the myosin from the network
                if sum((cond1+cond2).*(cond3+cond4)) + sum((cond5+cond6).*(cond7+cond8)) == 0
                    P.myosin_T2(aux_T3(lll)) = P.myosin_Tfree(aux_new(aux_c(aux_r)));
                    r_s(myosin2(aux_T3(lll)),:) = r_s_aux;
                    P.myosin_Tfree(aux_new(aux_c(aux_r))) = [];
                else
                    edges_s3 = [edges_s3;aux(aux2)];
                    aux_TT = [aux_TT;aux3];
                    aux_TT2 = [aux_TT2;aux4];
                end
            else
                disp(aux(aux2))
                edges_s3 = [edges_s3;aux(aux2)];
                aux_TT = [aux_TT;aux3];
                aux_TT2 = [aux_TT2;aux4];
            end
        end
        P.myosin_Tfree = [P.myosin_Tfree;P.myosin_T2(aux_TT2)];
        myosin2(aux_TT) = [];
        P.myosin_T2(aux_TT) = [];
    end
end