function [edges_s,edge_type,myosin2,P] = remove_small_myosin2(r_s,edges_s,edge_type,myosin2,P)
% check the length of the myosin rods
    aux = find(edge_type == 4);
    d = r_s(edges_s(aux,1),:) - r_s(edges_s(aux,2),:);
    d = sqrt(dot(d,d,2));
    aux_r = find(d < P.min_r2);
    %     remove those myosin rods that are shorter than the minimum length
    if ~isempty(aux_r)
        aux_rem = [];
        for l = 1:length(aux_r)
            aux_r1 = find(myosin2 == edges_s(aux(aux_r(l)),1));
            if ~isempty(aux_r1)
                aux_r2 = find(myosin2 == edges_s(aux(aux_r(l)),2));
            else
                aux_r1 = find(myosin2 == edges_s(aux(aux_r(l)),2));
                aux_r2 = find(myosin2 == edges_s(aux(aux_r(l)),1));
            end
            aux_rem = [aux_rem;aux_r1;aux_r2];
        end
        myosin2(aux_rem) = [];
        P.myosin_Tfree = [P.myosin_Tfree;P.myosin_T2(aux_rem)];
        P.myosin_T2(aux_rem) = [];
        edges_s(aux(aux_r),:) = [];
        edge_type(aux(aux_r)) = [];
    end
end