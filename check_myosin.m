function [r_s,myosin,edges_s2,P] = check_myosin(T_s,r_s,edges_s,edge_type,myosin,edges_s1,edges_s2,P)
% check if myosin linkers can attach to other spectrin triangles

% check the corresponding triangles
    aux_T = [find(T_s(P.myosin_T,1) == edges_s(edges_s1,2));
        find(T_s(P.myosin_T,2) == edges_s(edges_s1,2));
        find(T_s(P.myosin_T,3) == edges_s(edges_s1,2))];
          
    if ~isempty(aux_T)
        aux_TT = [];
        for lll = 1:length(aux_T)
%             if the edge corresponded to the triangle, then check whether
%             the other end of the myosin can attach to another triangle
%             within the allowed lengths of myosin 
            if ~isempty(intersect(edges_s(edges_s1,1),[T_s(P.myosin_T(aux_T(lll)),1);
                T_s(P.myosin_T(aux_T(lll)),2);T_s(P.myosin_T(aux_T(lll)),3)]))
            
                center_Tfree = (r_s(T_s(P.myosin_Tfree,1),:) + r_s(T_s(P.myosin_Tfree,2),:) +...
                    r_s(T_s(P.myosin_Tfree,3),:))./3;
                auxb = find(edge_type == 3);
                auxb3 = find(P.stress == edges_s(auxb(edges_s(auxb,2) == myosin(aux_T(lll))),1));

                aux_d = center_Tfree - r_s(P.stress(auxb3),:);
                aux_d = sqrt(dot(aux_d,aux_d,2));
                aux_new = find(aux_d < P.max_r & aux_d > P.min_r);

                aux = find(edge_type == 4);
  
                if ~isempty(aux_new)
%                     randomly select one of the traingles, trace the
%                     myosin linker and check if it intersects with the
%                     other myosins 

                    aux_c = 1:length(aux_new);
                    aux_r = randi(length(aux_c));

                    A = (r_s(edges_s([aux;auxb],2),2) - r_s(edges_s([aux;auxb],1),2))./...
                        (r_s(edges_s([aux;auxb],2),1) - r_s(edges_s([aux;auxb],1),1));
                    B = -A.*r_s(edges_s([aux;auxb],1),1) + r_s(edges_s([aux;auxb],1),2);
                    r_s_auxb = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
                    C = (r_s(P.stress(auxb3),2) - r_s_auxb(2))./...
                        (r_s(P.stress(auxb3),1) - r_s_auxb(1));
                    D = -C*r_s_auxb(1) + r_s_auxb(2);
                    Xb = (B-D)./(C-A);
                    Yb = (B.*C - D.*A)./(C-A);
                    aux_a = [find(A==inf);find(A==-inf)];
                    if ~isempty(aux_a)
                        aux_aa = r_s(edges_s([aux;auxb],1),1);
                        Xb(aux_a) = aux_aa(aux_a);
                        Yb(aux_a) = C*Xb(aux_a) + D;
                    end
                    if C == inf || C == -inf
                        Xb = r_s(P.stress(auxb3),1)*ones(size(B));
                        Yb = A.*Xb + B;
                    end

                    cond1b = (r_s(edges_s([aux;auxb],1),1)-P.d_inter < Xb).*(r_s(edges_s([aux;auxb],2),1)+P.d_inter > Xb);
                    cond2b = (r_s(edges_s([aux;auxb],2),1)-P.d_inter < Xb).*(r_s(edges_s([aux;auxb],1),1)+P.d_inter > Xb);
                    cond3b = (r_s_auxb(1)-P.d_inter < Xb).*(r_s(P.stress(auxb3),1) > Xb);
                    cond4b = (r_s_auxb(1)+P.d_inter > Xb).*(r_s(P.stress(auxb3),1) < Xb);

                    cond5b = (r_s(edges_s([aux;auxb],1),2)-P.d_inter < Yb).*(r_s(edges_s([aux;auxb],2),2)+P.d_inter > Yb);
                    cond6b = (r_s(edges_s([aux;auxb],2),2)-P.d_inter < Yb).*(r_s(edges_s([aux;auxb],1),2)+P.d_inter > Yb);
                    cond7b = (r_s_auxb(2)-P.d_inter < Yb).*(r_s(P.stress(auxb3),2) > Yb);
                    cond8b = (r_s_auxb(2)+P.d_inter > Yb).*(r_s(P.stress(auxb3),2) < Yb);
                    while (sum((cond1b+cond2b).*(cond3b+cond4b)) + sum((cond5b+cond6b).*(cond7b+cond8b)) > 0)...
                        && length(aux_c)>1 

                        aux_c(aux_r) = [];
                        aux_r = randi(length(aux_c));

                        r_s_auxb = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
                        C = (r_s(P.stress(auxb3),2) - r_s_auxb(2))./...
                        (r_s(P.stress(auxb3),1) - r_s_auxb(1));
                        D = -C*r_s_auxb(1) + r_s_auxb(2);
                        Xb = (B-D)./(C-A);
                        Yb = (B.*C - D.*A)./(C-A);

                        aux_a = [find(A==inf);find(A==-inf)];
                        if ~isempty(aux_a)
                            aux_aa = r_s(edges_s([aux;auxb],1),1);
                            Xb(aux_a) = aux_aa(aux_a);
                            Yb(aux_a) = C*Xb(aux_a) + D;
                        end
                        if C == inf || C == -inf
                            Xb = r_s(P.stress(auxb3),1)*ones(size(B));
                            Yb = A.*Xb + B;
                        end
                                
                        cond1b = (r_s(edges_s([aux;auxb],1),1)-P.d_inter < Xb).*(r_s(edges_s([aux;auxb],2),1)+P.d_inter > Xb);
                        cond2b = (r_s(edges_s([aux;auxb],2),1)-P.d_inter < Xb).*(r_s(edges_s([aux;auxb],1),1)+P.d_inter > Xb);
                        cond3b = (r_s_auxb(1)-P.d_inter < Xb).*(r_s(P.stress(auxb3),1) > Xb);
                        cond4b = (r_s_auxb(1)+P.d_inter > Xb).*(r_s(P.stress(auxb3),1) < Xb);

                        cond5b = (r_s(edges_s([aux;auxb],1),2)-P.d_inter < Yb).*(r_s(edges_s([aux;auxb],2),2)+P.d_inter > Yb);
                        cond6b = (r_s(edges_s([aux;auxb],2),2)-P.d_inter < Yb).*(r_s(edges_s([aux;auxb],1),2)+P.d_inter > Yb);
                        cond7b = (r_s_auxb(2)-P.d_inter < Yb).*(r_s(P.stress(auxb3),2) > Yb);
                        cond8b = (r_s_auxb(2)+P.d_inter > Yb).*(r_s(P.stress(auxb3),2) < Yb);
                    end
%                     check if the myosin linker does not intersect the other myosins and update its position, otherwise 
%                 remove the myosin from the network
                    if sum((cond1b+cond2b).*(cond3b+cond4b)) + sum((cond5b+cond6b).*(cond7b+cond8b)) == 0
                    
                        P.myosin_T(aux_T(lll)) = P.myosin_Tfree(aux_new(aux_c(aux_r)));
                        P.myosin_Tfree(aux_new(aux_c(aux_r))) = [];
                        r_s(myosin(aux_T(lll)),:) = r_s_auxb;
                    else
                        aux_TT = [aux_TT;aux_T(lll)];
                        edges_s2 = [edges_s2;auxb(aux_T(lll))];
                    end
                else
                    aux_TT = [aux_TT;aux_T(lll)];
                    edges_s2 = [edges_s2;auxb(aux_T(lll))];
                end
            end
        end
        P.myosin_Tfree = [P.myosin_Tfree;P.myosin_T(aux_TT)];
        myosin(aux_TT) = [];
        P.myosin_T(aux_TT) = [];
    end
end