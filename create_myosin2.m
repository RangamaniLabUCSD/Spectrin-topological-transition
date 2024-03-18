function [r_s,edges_s,edge_type,myosin2,P] = create_myosin2(r_s,T_s,edges_s,edge_type,myosin2,P)
% randomly choose a free spectrin triangle in the mesh 
    aux_r1 = randi(length(P.myosin_Tfree));
%     check the distance to the other free triangle centers
    center_Tfree = (r_s(T_s(P.myosin_Tfree,1),:) + r_s(T_s(P.myosin_Tfree,2),:) +...
        r_s(T_s(P.myosin_Tfree,3),:))./3;
    aux_d = center_Tfree - center_Tfree(aux_r1,:);
    aux_d = sqrt(dot(aux_d,aux_d,2));
%     select those triangles that are within the myosin minimum and maximum
%     lenght 
    aux_new = find(aux_d < P.max_r2 & aux_d > P.min_r2);
    
    if ~isempty(aux_new)
%       randomly select one of the selected triangles   
        aux_c = 1:length(aux_new);
        aux_r = randi(length(aux_c));
%         check if the edge expanding between the selected spectrin
%         triangles intersect other myosin rod (edge_type = 4) or myosin
%         linker (edge_type = 3)
        aux = find(edge_type == 4);

        auxb = find(edge_type == 3);
        
        if ~isempty(aux) || ~isempty(auxb)

            a = (r_s(edges_s([aux;auxb],2),2) - r_s(edges_s([aux;auxb],1),2))./...
                (r_s(edges_s([aux;auxb],2),1) - r_s(edges_s([aux;auxb],1),1));
            b = -a.*r_s(edges_s([aux;auxb],1),1) + r_s(edges_s([aux;auxb],1),2);
            r_s_aux = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
            c = (center_Tfree(aux_r1,2) - r_s_aux(2))./...
                (center_Tfree(aux_r1,1) - r_s_aux(1));
            d = -c*r_s_aux(1) + r_s_aux(2);
            X = (b-d)./(c-a);
            Y = (b.*c - d.*a)./(c-a);
%             correct for vertical lines
            aux_a = [find(a==inf);find(a==-inf)];
            if ~isempty(aux_a)
                aux_aa = r_s(edges_s([aux;auxb],1),1);
                X(aux_a) = aux_aa(aux_a);
                Y(aux_a) = c*X(aux_a) + d;
            end
            if c == inf || c == -inf
                X = center_Tfree(aux_r1,1)*ones(size(b));
                Y = a.*X + b;
            end
%             check that the myosins intersect in their length +- P.d_inter
            cond1 = (r_s(edges_s([aux;auxb],1),1)-P.d_inter < X).*(r_s(edges_s([aux;auxb],2),1)+P.d_inter > X);
            cond2 = (r_s(edges_s([aux;auxb],2),1)-P.d_inter < X).*(r_s(edges_s([aux;auxb],1),1)+P.d_inter > X);
            cond3 = (r_s_aux(1)-P.d_inter < X).*(center_Tfree(aux_r1,1)+P.d_inter > X);
            cond4 = (r_s_aux(1)+P.d_inter > X).*(center_Tfree(aux_r1,1)-P.d_inter < X);
            cond5 = (r_s(edges_s([aux;auxb],1),2)-P.d_inter < Y).*(r_s(edges_s([aux;auxb],2),2)+P.d_inter > Y);
            cond6 = (r_s(edges_s([aux;auxb],2),2)-P.d_inter < Y).*(r_s(edges_s([aux;auxb],1),2)+P.d_inter > Y);
            cond7 = (r_s_aux(2)-P.d_inter < Y).*(center_Tfree(aux_r1,2)+P.d_inter > Y);
            cond8 = (r_s_aux(2)+P.d_inter > Y).*(center_Tfree(aux_r1,2)-P.d_inter < Y);
%          if they intersect, take another free triangle  
            while (sum((cond1+cond2).*(cond3+cond4)) + sum((cond5+cond6).*(cond7+cond8)) > 0)...
                                && length(aux_c)>1 
              
                aux_c(aux_r) = [];
                aux_r = randi(length(aux_c));

                r_s_aux = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
                c = (center_Tfree(aux_r1,2) - r_s_aux(2))./...
                (center_Tfree(aux_r1,1) - r_s_aux(1));
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
                    X = center_Tfree(aux_r1,1)*ones(size(b));
                    Y = a.*X + b;
                end
                cond1 = (r_s(edges_s([aux;auxb],1),1)-P.d_inter < X).*(r_s(edges_s([aux;auxb],2),1)+P.d_inter > X);
                cond2 = (r_s(edges_s([aux;auxb],2),1)-P.d_inter < X).*(r_s(edges_s([aux;auxb],1),1)+P.d_inter > X);
                cond3 = (r_s_aux(1)-P.d_inter < X).*(center_Tfree(aux_r1,1)+P.d_inter > X);
                cond4 = (r_s_aux(1)+P.d_inter > X).*(center_Tfree(aux_r1,1)-P.d_inter < X);
                cond5 = (r_s(edges_s([aux;auxb],1),2)-P.d_inter < Y).*(r_s(edges_s([aux;auxb],2),2)+P.d_inter > Y);
                cond6 = (r_s(edges_s([aux;auxb],2),2)-P.d_inter < Y).*(r_s(edges_s([aux;auxb],1),2)+P.d_inter > Y);
                cond7 = (r_s_aux(2)-P.d_inter < Y).*(center_Tfree(aux_r1,2)+P.d_inter > Y);
                cond8 = (r_s_aux(2)+P.d_inter > Y).*(center_Tfree(aux_r1,2)-P.d_inter < Y);
            end
%  check whether there was a myosin that does not intersect others and add
%  it to the list 
            if sum((cond1+cond2).*(cond3+cond4)) + sum((cond5+cond6).*(cond7+cond8)) == 0
                P.myosin_T2 =  [P.myosin_T2;P.myosin_Tfree(aux_r1);P.myosin_Tfree(aux_new(aux_c(aux_r)))];
                P.myosin_Tfree([aux_new(aux_c(aux_r));aux_r1]) = [];
                edges_s = [edges_s;size(r_s,1)+(1:2)];
                edge_type = [edge_type;4];
                myosin2 = [myosin2;size(r_s,1)+(1:2)'];
                r_s = [r_s;sum(r_s(T_s(P.myosin_T2(end-1),:),:))./3;sum(r_s(T_s(P.myosin_T2(end),:),:))./3];
            end
        else
%  if there are no other myosins, just add the new one             
            aux_r = randi(length(aux_new));
            P.myosin_T2 =  [P.myosin_T2;P.myosin_Tfree(aux_r1);P.myosin_Tfree(aux_new(aux_r))];
            P.myosin_Tfree([aux_new(aux_r);aux_r1]) = [];
            edges_s = [edges_s;size(r_s,1)+(1:2)];
            edge_type = [edge_type;4];
            myosin2 = [myosin2;size(r_s,1)+(1:2)'];
            r_s = [r_s;sum(r_s(T_s(P.myosin_T2(end-1),:),:))./3;sum(r_s(T_s(P.myosin_T2(end),:),:))./3];
        end
    end
end