function [edges_s,edge_type,myosin,P] = remove_small_myosin(r_s,edges_s,edge_type,myosin,P)   
% check the length of the myosin linkers 
    aux = find(edge_type == 3);
    d = r_s(edges_s(aux,1),:) - r_s(edges_s(aux,2),:);
    d = sqrt(dot(d,d,2));
    aux_r = find(d < P.min_r);
%     remove those myosin linkers that are shorter than the minimum length
    if ~isempty(aux_r)
        aux_rem1 = [];
        for l = 1:length(aux_r)
            aux_r1 = find(myosin == edges_s(aux(aux_r(l)),2));
            aux_rem1 = [aux_rem1;aux_r1];
        end
        myosin(aux_rem1) = [];
        P.myosin_Tfree = [P.myosin_Tfree;P.myosin_T(aux_rem1)];
        P.myosin_T(aux_rem1) = [];
        edges_s(aux(aux_r),:) = [];
        edge_type(aux(aux_r)) = [];
    end
end