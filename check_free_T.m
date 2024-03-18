function P = check_free_T(T_s,edges_s,edge_type,P)
% check which triangles formed by spectrin bundles (edge_type = 0 ) do not have myosin
% attached.  
    P.myosin_Tfree = [];
    for l = 1:size(T_s,1)
        aux_e1 = find(edges_s(edge_type == 0,1) == T_s(l,1) & edges_s(edge_type == 0,2) == T_s(l,2));
        if ~isempty(aux_e1)
            aux_e2 = find(edges_s(edge_type == 0,1) == T_s(l,2) & edges_s(edge_type == 0,2) == T_s(l,3));
            if ~isempty(aux_e2)
                aux_e3 = find(edges_s(edge_type == 0,1) == T_s(l,3) & edges_s(edge_type == 0,2) == T_s(l,1));
                if ~isempty(aux_e3)
                    P.myosin_Tfree = [P.myosin_Tfree;l];
                else
                    aux_e3 = find(edges_s(edge_type == 0,2) == T_s(l,3) & edges_s(edge_type == 0,1) == T_s(l,1));
                    if ~isempty(aux_e3)
                        P.myosin_Tfree = [P.myosin_Tfree;l];
                    end
                end
            else
                aux_e2 = find(edges_s(edge_type == 0,2) == T_s(l,2) & edges_s(edge_type == 0,1) == T_s(l,3));
                if ~isempty(aux_e2)
                    aux_e3 = find(edges_s(edge_type == 0,1) == T_s(l,3) & edges_s(edge_type == 0,2) == T_s(l,1));
                    if ~isempty(aux_e3)
                        P.myosin_Tfree = [P.myosin_Tfree;l];
                    else
                        aux_e3 = find(edges_s(edge_type == 0,2) == T_s(l,3) & edges_s(edge_type == 0,1) == T_s(l,1));
                        if ~isempty(aux_e3)
                            P.myosin_Tfree = [P.myosin_Tfree;l];
                        end
                    end
                end
            end
        else
            aux_e1 = find(edges_s(edge_type == 0,2) == T_s(l,1) & edges_s(edge_type == 0,1) == T_s(l,2));
            if ~isempty(aux_e1)
                aux_e2 = find(edges_s(edge_type == 0,1) == T_s(l,2) & edges_s(edge_type == 0,2) == T_s(l,3));
                if ~isempty(aux_e2)
                    aux_e3 = find(edges_s(edge_type == 0,1) == T_s(l,3) & edges_s(edge_type == 0,2) == T_s(l,1));
                    if ~isempty(aux_e3)
                        P.myosin_Tfree = [P.myosin_Tfree;l];
                    else
                        aux_e3 = find(edges_s(edge_type == 0,2) == T_s(l,3) & edges_s(edge_type == 0,1) == T_s(l,1));
                        if ~isempty(aux_e3)
                            P.myosin_Tfree = [P.myosin_Tfree;l];
                        end
                    end
                else
                    aux_e2 = find(edges_s(edge_type == 0,2) == T_s(l,2) & edges_s(edge_type == 0,1) == T_s(l,3));
                    if ~isempty(aux_e2)
                        aux_e3 = find(edges_s(edge_type == 0,1) == T_s(l,3) & edges_s(edge_type == 0,2) == T_s(l,1));
                        if ~isempty(aux_e3)
                            P.myosin_Tfree = [P.myosin_Tfree;l];
                        else
                            aux_e3 = find(edges_s(edge_type == 0,2) == T_s(l,3) & edges_s(edge_type == 0,1) == T_s(l,1));
                            if ~isempty(aux_e3)
                                P.myosin_Tfree = [P.myosin_Tfree;l];
                            end
                        end
                    end
                end
            end
        end
    end
end